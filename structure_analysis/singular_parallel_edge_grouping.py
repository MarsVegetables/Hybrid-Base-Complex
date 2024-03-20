from tqdm import tqdm
from mesh_structure.mesh_common import reduce_id_groups, check_common_elements


class SingularParallelEdgeGrouping:
    def __init__(self, singularity_graph):
        self.sg = singularity_graph
        self.region_eids = self.sg.region_eids
        self.mesh = self.sg.mesh
        self.edge_singular_edge_map = self.sg.edge_singular_edge_map
        self.singular_es = self.sg.singular_es
        # matrix is represented by dict
        self.mesh_parallel_matrix = {}
        self.singular_edge_parallel_matrix = {}
        self.singular_edge_parallel_ratio_matrix = {}
        self.visited_eids = set()
        #
        self.singular_edge_group_id = {}

    def find_parallel_singular_edges(self, singular_edge_id):
        se = self.singular_es[singular_edge_id]
        for mesh_eid in se.es_link:
            self.find_all_parallel_mesh_edges(mesh_eid)

    def find_all_parallel_mesh_edges(self, mesh_eid):
        mesh_eid_stack = [mesh_eid]
        parallel_mesh_eids = set()
        grouped_singular_eids = set()
        while mesh_eid_stack:
            eid = mesh_eid_stack.pop()
            if eid in self.visited_eids:
                continue
            self.visited_eids.add(eid)
            parallel_mesh_eids.add(eid)
            # add edge's parallel edges
            edge = self.mesh.Edges[eid]
            region_parallel_edges = check_common_elements(edge.parallels, self.region_eids)
            mesh_eid_stack.extend(region_parallel_edges)
            # add singular edge into matrix
            if eid not in self.edge_singular_edge_map:
                continue
            se_ids = self.edge_singular_edge_map[eid]
            grouped_singular_eids.update(se_ids)
        # complete parallel matrix
        for eid in parallel_mesh_eids:
            if eid not in self.mesh_parallel_matrix:
                self.mesh_parallel_matrix[eid] = set()
            self.mesh_parallel_matrix[eid].update(parallel_mesh_eids)
            self.mesh_parallel_matrix[eid].remove(eid)
        # find all parallel singular edges
        for seid in grouped_singular_eids:
            if seid not in self.singular_edge_parallel_matrix:
                self.singular_edge_parallel_matrix[seid] = {}
            for seid_2 in grouped_singular_eids:
                if seid_2 == seid:
                    continue
                if seid_2 not in self.singular_edge_parallel_matrix[seid]:
                    self.singular_edge_parallel_matrix[seid][seid_2] = 0
                self.singular_edge_parallel_matrix[seid][seid_2] += 1

    # step 1
    def complete_parallel_matrix(self):
        self.visited_eids = set()
        self.mesh_parallel_matrix = {}
        self.singular_edge_parallel_matrix = {}
        for i, se in tqdm(enumerate(self.singular_es)):
            self.find_parallel_singular_edges(i)

    # step 2
    def compute_singular_edge_parallel_ratio(self):
        self.singular_edge_parallel_ratio_matrix = {sid: {} for sid in self.singular_edge_parallel_matrix}
        for seid in self.singular_edge_parallel_matrix:
            singular_edge_n = len(self.singular_es[seid].es_link)
            for seid_2 in self.singular_edge_parallel_matrix[seid]:
                parallel_edge_number = self.singular_edge_parallel_matrix[seid][seid_2]
                if seid_2 not in self.singular_edge_parallel_ratio_matrix[seid]:
                    self.singular_edge_parallel_ratio_matrix[seid][seid_2] = 0
                self.singular_edge_parallel_ratio_matrix[seid][seid_2] = parallel_edge_number / singular_edge_n

    # step 3
    def group_singular_edges(self, parallel_ratio_threshold):
        edge_groups = []
        for seid in self.singular_edge_parallel_ratio_matrix:
            g = [seid]
            row = self.singular_edge_parallel_ratio_matrix[seid]
            for seid_2 in row:
                ratio = row[seid_2]
                if ratio < parallel_ratio_threshold:
                    continue
                g.append(seid_2)
            edge_groups.append(g)
        cleaned_groups = reduce_id_groups(edge_groups)
        return cleaned_groups

    def grouping(self, parallel_ratio_threshold):
        self.complete_parallel_matrix()
        self.compute_singular_edge_parallel_ratio()
        groups = self.group_singular_edges(parallel_ratio_threshold)
        self.singular_edge_group_id = {}
        for i, g in enumerate(groups):
            for se_id in g:
                self.singular_edge_group_id[se_id] = i
