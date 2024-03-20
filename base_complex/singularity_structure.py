# https://github.dev/Cotrik/CotrikMesh/blob/master/libcotrik/src/Mesh.cpp
# https://github.com/Cotrik/CotrikMesh/blob/master/libcotrik/src/BaseComplex.cpp
from base_complex.base_complex_elements import SingularityVert, SingularityEdge


class SingularityStructure:
    def __init__(self, mesh, iteration_num=0):
        self.mesh = mesh
        self.init_attrs(iteration_num)
        self.build()

    def init_attrs(self, iteration_num=0):
        self.singular_Vs = []  # singular vert list
        self.singular_Es = []  # singular edge list
        # self.isComponent = is_component
        self.iter_num = iteration_num

    # ToDo : consider if two edge contain same valence and share cell
    def get_next_singular_edge_id(self, vid, eid):
        nextEdgeId = -1
        cur_edge = self.mesh.Edges[eid]
        num_1 = 0
        num_2 = 0
        for neid in self.mesh.Vertices[vid].neighboring_Eids:
            if neid == eid:
                continue
            next_edge = self.mesh.Edges[neid]
            # same valece
            if len(next_edge.neighboring_Cids) == len(cur_edge.neighboring_Cids):
                num_1 += 1
                nextEdgeId = neid
            elif next_edge.is_singular:
                num_2 += 1
        if (num_1 == 1) and (num_2 == 0):
            return nextEdgeId
        return -1

    def trace_edge(self, vid, eid, singular_Vids, singular_Eids):
        cur_eid = eid
        if not self.mesh.Edges[cur_eid].visited:
            singular_Eids.append(cur_eid)
            self.mesh.Edges[cur_eid].visited = True

        cur_vid = vid
        singular_Vids.append(cur_vid)
        next_eid = self.get_next_singular_edge_id(cur_vid, cur_eid)
        while next_eid != -1:
            singular_Eids.append(next_eid)
            cur_eid = next_eid
            self.mesh.Edges[cur_eid].visited = True
            if self.mesh.Edges[cur_eid].Vids[0] == cur_vid:
                cur_vid = self.mesh.Edges[cur_eid].Vids[1]
            else:
                cur_vid = self.mesh.Edges[cur_eid].Vids[0]
            singular_Vids.append(cur_vid)
            if cur_eid == eid:
                return True
            next_eid = self.get_next_singular_edge_id(cur_vid, cur_eid)
        return False

    def add_circular_singular_edge(self, leftVids, leftEids):
        singularEdgeId = len(self.singular_Es)
        se = SingularityEdge(singularEdgeId)
        se.start_end_vids = [-1, -1]
        for leftEid in leftEids:
            se.es_link.append(leftEid)
            self.mesh.Edges[leftEid].seid = singularEdgeId
        se.vs_link.extend(leftVids)
        se.is_boundary = self.mesh.Edges[leftEids[0]].is_boundary
        se.is_circle = True
        self.singular_Es.append(se)

    def add_singular_edge(self, leftVids, leftEids, rightVids, rightEids, left_svid, right_svid):
        new_seid = len(self.singular_Es)
        se = SingularityEdge(new_seid)
        se.start_end_vids = [left_svid, right_svid]
        se.es_link.extend(leftEids)
        se.es_link.extend(rightEids)
        for eid in se.es_link:
            self.mesh.Edges[eid].seid = new_seid
        se.vs_link.extend(leftVids)
        se.vs_link.extend(rightVids)
        se.is_boundary = self.mesh.Edges[se.es_link[0]].is_boundary
        self.singular_Es.append(se)
        return True

    def add_singular_vert(self, vid):
        for sv in self.singular_Vs:
            if sv.mesh_vid == vid:
                return sv.element_id
        svid = len(self.singular_Vs)
        sv = SingularityVert(svid, vid)
        sv.is_boundary = self.mesh.Vertices[vid].is_boundary
        self.singular_Vs.append(sv)
        return svid

    def extract_singular_v_and_e(self):
        print("extracting singular vertices and edges ...")
        self.mesh.clean_visited()
        self.mesh.mark_singularity_elements()
        for eid, edge in self.mesh.Edges.items():
            if edge.visited or not edge.is_singular:
                continue
            # =============================================================
            # left
            # =============================================================
            leftVids = []
            leftEids = []
            isCircle = self.trace_edge(edge.Vids[0], eid, leftVids, leftEids)
            # =============================================================
            # circle
            # if left direction is a circle
            # =============================================================
            if isCircle:
                self.add_circular_singular_edge(leftVids, leftEids)
                continue
            leftSingularVid = self.add_singular_vert(leftVids[-1])
            # =============================================================
            # right
            # =============================================================
            rightVids = []
            rightEids = []
            isCircle = self.trace_edge(edge.Vids[1], eid, rightVids, rightEids)
            rightSingularVid = self.add_singular_vert(rightVids[-1])
            self.add_singular_edge(leftVids, leftEids, rightVids, rightEids, leftSingularVid, rightSingularVid)
        print("done ...")

    def build_singularity_connectivity(self):
        print("building singularity vertex connextivities ...")
        for se in self.singular_Es:
            v0 = se.start_end_vids[0]
            v1 = se.start_end_vids[1]
            # not circular
            if v0 != v1:
                self.singular_Vs[v0].neighboring_Vids.append(v1)
                self.singular_Vs[v0].neighboring_Eids.append(se.element_id)
                self.singular_Vs[v1].neighboring_Vids.append(v0)
                self.singular_Vs[v1].neighboring_Eids.append(se.element_id)
        print("done ...")

    def build(self):
        self.extract_singular_v_and_e()
        self.build_singularity_connectivity()

    def clean_singular_element_visited(self):
        for se in self.singular_Es:
            se.visited = False
        for sv in self.singular_Vs:
            sv.visited = False

    def get_vertex_list(self):
        vertexList = []
        for sv in self.singular_Vs:
            coord = self.mesh.Vertices[sv.mesh_id].xyz()
            vertexList.append(coord)
        return vertexList


