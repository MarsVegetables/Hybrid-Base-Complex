from mesh_structure.mesh_common import check_common_elements, reduce_id_groups
from tqdm import tqdm


# matching chain edge grouping
class MaximumMatchingFinder:
    def __init__(self, mesh, edge_scope, cell_scope):
        self.mesh = mesh
        # parallel edges
        self.parallel_edges = edge_scope
        self.cell_scope = cell_scope
        self.one_ring_cell_scope = set()
        # include all cells neighboring with parallel edges
        self.complete_edge_one_ring_cells()
        # sheet is built by connect hex cells
        # intersecting is also caused by hex
        # self.one_ring_cells_include_non_hex = set()
        # region edge
        self.region_edges = set()
        self.complete_region_edges()
        self.visited_eids = set()
        self.visited_cids = set()

    def complete_edge_one_ring_cells(self):
        self.one_ring_cell_scope = set()
        for eid in self.parallel_edges:
            edge = self.mesh.Edges[eid]
            ncids = edge.neighboring_Cids
            self.one_ring_cell_scope.update(ncids)

    def complete_region_edges(self):
        self.region_edges = set()
        for cid in self.cell_scope:
            cell = self.mesh.Cells[cid]
            self.region_edges.update(cell.Eids)

    def decompose_self_intersecting_sheet_matching(self):
        self.visited_eids = set()
        all_matching_edge_sets = []
        all_matching_cell_sets = []
        for eid in self.parallel_edges:
            if eid in self.visited_eids:
                continue
            ms, mc = self.find_maximum_matching(eid)
            if ms:
                all_matching_edge_sets.append(ms)
                all_matching_cell_sets.append(mc)
        self.visited_eids = set()
        return all_matching_edge_sets, all_matching_cell_sets

    # maximum matching
    # edges in the set does not share any vertices
    def find_maximum_matching(self, eid):
        visited_vids = set()
        matching_set = set()
        candidate_edges = [eid]
        while candidate_edges:
            cur_eid = candidate_edges.pop()
            if cur_eid in self.visited_eids:
                continue
            edge = self.mesh.Edges[cur_eid]
            if check_common_elements(edge.Vids, visited_vids):
                continue
            visited_vids.update(edge.Vids)
            self.visited_eids.add(cur_eid)
            matching_set.add(cur_eid)
            next_parallel_edges = check_common_elements(self.parallel_edges, edge.parallels)
            candidate_edges.extend(next_parallel_edges)
        cell_set = set()
        for eid in matching_set:
            edge = self.mesh.Edges[eid]
            cids = check_common_elements(edge.neighboring_Cids, self.cell_scope)
            cell_set.update(cids)
        return matching_set, cell_set

    def edge_visited(self, eid):
        edge = self.mesh.Edges[eid]
        pe = check_common_elements(edge.parallels, self.parallel_edges)
        visited_pe = check_common_elements(pe, self.visited_eids)
        if visited_pe == pe:
            return True
        return False

    def update_layer_parallel_edges(self, parallel_eids, layer_cells):
        res = set()
        edge_scope = set()
        for cid in layer_cells:
            cell = self.mesh.Cells[cid]
            edge_scope.update(cell.Eids)
        valid_parallel_eids = check_common_elements(parallel_eids, edge_scope)
        eid_stack = list(valid_parallel_eids)
        while eid_stack:
            cur_eid = eid_stack.pop()
            if cur_eid in res:
                continue
            if cur_eid not in edge_scope:
                continue
            edge = self.mesh.Edges[cur_eid]
            valid_e = check_common_elements(valid_parallel_eids, edge.parallels)
            eid_stack.extend(valid_e)
            remain_e = set(edge.parallels) - valid_e
            for pe in remain_e:
                if self.is_parallel_edge_valid(pe, valid_parallel_eids, layer_cells):
                    eid_stack.append(pe)
            res.add(cur_eid)
        return res

    # layer contain self-parallel sheets
    # decompose based on common cell condition
    def decompose_self_intersecting_sheet_layer(self):
        print("decomposing self-intersecting sheet ...")
        # self.visited_cids = set()
        self.visited_eids = set()
        all_layer_edge_sets = []
        all_layer_cell_sets = []
        for eid in tqdm(self.parallel_edges):
            if self.edge_visited(eid):
                continue
            ms, layer_cids = self.find_maximum_layer(eid)
            if not ms:
                continue
            layer_parallel_eids = self.update_layer_parallel_edges(ms, layer_cids)
            self.visited_eids.update(layer_parallel_eids)
            all_layer_edge_sets.append(layer_parallel_eids)
            all_layer_cell_sets.append(layer_cids)
        print("decomposing self-intersecting sheet ... done")
        return all_layer_edge_sets, all_layer_cell_sets

    # parallel edge
    # neighoring with a cell in the layer
    # and does not make intersecting with existing parallel edges
    def is_parallel_edge_valid(self, parallel_eid, layer_parallel_eids, layer_cids):
        p_edge = self.mesh.Edges[parallel_eid]
        if layer_cids and (not check_common_elements(p_edge.neighboring_Cids, layer_cids)):
            return False
        in_set_eids = check_common_elements(p_edge.neighboring_Eids, layer_parallel_eids)
        # if neighboring edge is intersecting
        for eid in in_set_eids:
            edge = self.mesh.Edges[eid]
            common_cids = check_common_elements(p_edge.neighboring_Cids, edge.neighboring_Cids)
            common_cids = check_common_elements(common_cids, self.cell_scope)
            if common_cids:
                return False
        return True

    def find_maximum_layer(self, eid):
        layer_set = set()
        layer_cells = set()
        candidate_edges = [eid]
        while candidate_edges:
            cur_eid = candidate_edges.pop()
            if cur_eid in layer_set:
                continue
            if self.edge_visited(cur_eid):
                continue
            edge = self.mesh.Edges[cur_eid]
            # not intersecting with current parallel edges
            if not self.is_parallel_edge_valid(cur_eid, layer_set, layer_cells):
                layer_cells = layer_cells - set(edge.neighboring_Cids)
                continue
            # for valid edge, add its neighbors into the layer
            edge_available_ncids = check_common_elements(edge.neighboring_Cids, self.cell_scope)
            # valid next neighbor
            if layer_cells:
                valid_neighbor_cells = set()
                valid_edge_set = set()
                for ncid in edge_available_ncids:
                    ncell = self.mesh.Cells[ncid]
                    if check_common_elements(ncell.neighboring_Cids, layer_cells):
                        valid_neighbor_cells.add(ncid)
                        valid_edge_set.update(ncell.Eids)
            else:
                valid_neighbor_cells = edge_available_ncids
                valid_edge_set = edge.parallels
            layer_cells.update(valid_neighbor_cells)
            # try next parallel edges
            edge_pe = check_common_elements(self.parallel_edges, edge.parallels)
            edge_pe = check_common_elements(edge_pe, valid_edge_set)
            candidate_edges.extend(edge_pe)
            # add the edge in the layer edge set
            layer_set.add(cur_eid)
        return layer_set, layer_cells
