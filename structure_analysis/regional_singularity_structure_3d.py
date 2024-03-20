# extract singular graph of a self-intersect sheet self-parallel sheet
# perfect sheet and type 1 sheet does not have complex volume configuration
from mesh_structure.mesh_common import check_common_elements
from base_complex.base_complex_elements import SingularityEdge
import random


# build singularity graph based on singular edges
# trace singular edge and assign color to each singular edge
class SingularityStructure3D:
    def __init__(self, mesh, in_scope_cells=None):
        self.mesh = mesh
        # elements
        self.cell_scope = set()
        # fids in scope need to determine the boundary elements
        self.in_scope_fids = set()
        self.in_scope_eids = set()
        self.in_scope_vids = set()

        self.boundary_fids = set()
        self.boundary_eids = set()
        self.boundary_vids = set()

        self.color_num = 0

        self.visited_eids = set()

        # just element ids
        self.singular_vs = set()
        # singular edge object list
        self.singular_es = []

        self.vertex_singular_edge_map = {}

        # avoid edge set
        self.avoid_edge_set = set()

        self.update_in_scope_elements(in_scope_cells)

    # avoid edge set is parallel edge set
    def update_avoid_edge_set(self, eids):
        self.avoid_edge_set.update(eids)

    def update_in_scope_elements(self, cids):
        if cids:
            self.cell_scope = cids
        else:
            self.cell_scope = self.mesh.Cells
        self.update_in_scope_fids()
        self.update_in_scope_eids()
        self.update_in_scope_vids()
        self.find_boundary_elements()

    def update_in_scope_fids(self):
        for cid in self.cell_scope:
            cell = self.mesh.Cells[cid]
            self.in_scope_fids.update(cell.Fids)

    def update_in_scope_eids(self):
        for cid in self.cell_scope:
            cell = self.mesh.Cells[cid]
            self.in_scope_eids.update(cell.Eids)

    def update_in_scope_vids(self):
        for cid in self.cell_scope:
            cell = self.mesh.Cells[cid]
            self.in_scope_vids.update(cell.Vids)

    def find_boundary_elements(self):
        for fid in self.in_scope_fids:
            face = self.mesh.Faces[fid]
            ncids = check_common_elements(face.neighboring_Cids, self.cell_scope)
            if len(ncids) == 1:
                continue
            self.boundary_fids.add(fid)
            self.boundary_eids.update(face.Eids)
            self.boundary_vids.update(face.Vids)

    def build_singularity_graph(self):
        self.extract_singular_edges()
        self.extract_singular_vertices()
        self.build_singular_edge_neighbor_relation()

    def add_singular_edges(self, vert_link, edge_link):
        if len(edge_link) == 0:
            return False
        new_seid = len(self.singular_es)
        se = SingularityEdge(new_seid)
        if vert_link[0] == vert_link[-1]:
            se.start_end_vids = []
            se.is_circle = True
        else:
            se.start_end_vids = [vert_link[0], vert_link[-1]]
            se.is_circle = False
        se.es_link = edge_link
        se.vs_link = vert_link
        se.is_boundary = True if se.es_link[0] in self.boundary_eids else False
        self.singular_es.append(se)

        # build vert and singular edge map
        for v in se.vs_link:
            if v in self.vertex_singular_edge_map:
                self.vertex_singular_edge_map[v].add(new_seid)
            else:
                self.vertex_singular_edge_map[v] = {new_seid}
        return True

    def extract_singular_edges(self):
        print("extracting 3D volume singular edges ...")
        seids = set()
        for eid in self.in_scope_eids:
            edge = self.mesh.Edges[eid]
            is_singular = False
            cids = check_common_elements(edge.neighboring_Cids, self.cell_scope)
            if len(cids) != 4 and (edge not in self.boundary_eids):
                is_singular = True
            if len(cids) != 2 and (edge in self.boundary_eids):
                is_singular = True
            if not is_singular:
                continue
            seids.add(eid)
        self.visited_eids = set()
        for eid in seids:
            if eid in self.avoid_edge_set:
                continue
            edge = self.mesh.Edges[eid]
            # left link
            # vx <-- xxxx -- v0 <-- eid -- v1
            left_vids, left_eids = self.trace_edge(edge.Vids[1], eid)
            if left_vids[0] == left_vids[-1]:
                self.add_singular_edges(left_vids, left_eids)
                continue
            full_vids = left_vids
            full_eids = left_eids
            # right link
            # v0 -- eid --> v1 -- next --> vx
            self.visited_eids.remove(eid)
            right_vids, right_eids = self.trace_edge(edge.Vids[0], eid)
            if len(right_eids) > 0:
                left_vids.reverse()
                left_eids.reverse()
                full_vids = left_vids[:-2]
                full_eids = left_eids[:-1]
                full_vids.extend(right_vids)
                full_eids.extend(right_eids)
            self.add_singular_edges(full_vids, full_eids)
        print("done ...")

    # in 3d space, if the edge contain more than one next edge
    # the edge is singularity
    def trace_next_edge(self, vid, eid):
        # next element
        next_eid = self.mesh.get_next_edge(vid, eid)
        if len(next_eid) > 1:
            return []
        return next_eid

    # give an edge, trace to both direction
    def trace_edge(self, vid, eid):
        # sequence of edge is matter
        edge_link = []
        vert_link = [vid]
        next_eid = [eid]
        next_vid = vid
        while next_eid:
            cur_eid = next_eid.pop()
            if cur_eid in edge_link:
                continue
            if cur_eid in self.visited_eids:
                continue
            if cur_eid not in self.in_scope_eids:
                continue
            edge_link.append(cur_eid)

            self.visited_eids.add(cur_eid)

            edge = self.mesh.Edges[cur_eid]
            for evid in edge.Vids:
                if evid != next_vid:
                    next_vid = evid
                    break
            vert_link.append(next_vid)
            next_eid = self.trace_next_edge(next_vid, cur_eid)
        return vert_link, edge_link

    def add_singular_vertices(self, vids):
        self.singular_vs.update(vids)

    def extract_singular_vertices(self):
        print("extracting 3D volume singular vertices ...")
        # find high valence vertices
        svids = set()
        for se in self.singular_es:
            # is singularity vert
            svids.update(se.start_end_vids)
        self.add_singular_vertices(svids)
        print("done ...")

    def build_singular_edge_neighbor_relation(self):
        for vid in self.vertex_singular_edge_map:
            nsids = self.vertex_singular_edge_map[vid]
            for se_id in nsids:
                tmp = set(self.singular_es[se_id].neighboring_Eids)
                tmp.update(nsids)
                tmp.remove(se_id)
                self.singular_es[se_id].neighboring_Eids = tmp
        self.determine_color_number()
        self.assign_singular_edge_color()

    def determine_color_number(self):
        color_number = 0
        for se in self.singular_es:
            if len(se.neighboring_Eids) + 2 > color_number:
                color_number = len(se.neighboring_Eids) + 2
        self.color_num = color_number
        return color_number

    def assign_singular_edge_color(self):
        for se in self.singular_es:
            neighbor_colors = set()
            for n_se_id in se.neighboring_Eids:
                cc = self.singular_es[n_se_id]
                neighbor_colors.add(cc.color)
            candidate_colors = set()
            for i in range(self.color_num):
                # keep rad color reserved
                if i == 0:
                    continue
                if i not in neighbor_colors:
                    candidate_colors.add(i)
            picked_color = random.choice(list(candidate_colors))
            se.color = picked_color
