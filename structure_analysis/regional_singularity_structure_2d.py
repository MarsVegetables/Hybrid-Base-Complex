from mesh_structure.mesh_common import check_common_elements
from base_complex.base_complex_elements import SingularityEdge
import random
from copy import deepcopy


# analysis structure of perfect sheets and type-1 imperfect sheets
class SingularityStructure2D:
    def __init__(self, mesh, in_scope_cells=None, in_scope_faces=None, avoid_eids=None):
        self.mesh = mesh

        # cids
        self.input_cell_scope = set()
        self.extend_cells = set()
        self.cell_scope = set()

        # fids
        self.input_face_scope = set()
        self.face_scope = set()
        # extend face scope contain faces from non hex element
        self.extend_faces = set()

        # usually are parallel edges in a sheet
        self.avoid_edge_set = set()
        self.update_avoid_edge_set(avoid_eids)

        self.in_scope_vids = set()
        self.in_scope_eids = set()

        self.boundary_vids = set()
        self.boundary_eids = set()

        self.color_num = 0

        self.update_cell_scope(in_scope_cells)
        self.update_face_scope(in_scope_faces)

        self.visited_eids = set()

        # just element ids
        self.singular_vs = set()
        self.redundant_singular_vs = set()
        # singular edge object list
        self.singular_es = []

        self.vertex_singular_edge_map = {}

    # avoid edge set is parallel edge set
    def update_avoid_edge_set(self, eids):
        if not eids:
            return
        self.avoid_edge_set.update(eids)

    def find_boundary_edges_and_vertices(self):
        for eid in self.in_scope_eids:
            edge = self.mesh.Edges[eid]
            nfids = check_common_elements(edge.neighboring_Fids, self.face_scope)
            if len(nfids) != 2:
                self.boundary_eids.add(eid)
                self.boundary_vids.update(edge.Vids)

    def update_face_scope(self, fids):
        if fids:
            self.input_face_scope = deepcopy(fids)
            self.face_scope = fids
        else:
            self.face_scope = self.mesh.Faces
        self.extend_one_ring_faces()
        self.update_in_scope_vids()
        self.update_in_scope_eids()
        self.find_boundary_edges_and_vertices()

    def update_in_scope_vids(self):
        for fid in self.face_scope:
            face = self.mesh.Faces[fid]
            self.in_scope_vids.update(face.Vids)

    def update_in_scope_eids(self):
        for fid in self.face_scope:
            face = self.mesh.Faces[fid]
            self.in_scope_eids.update(face.Eids)

    def build_singularity_graph(self):
        self.extract_singular_vertices()
        self.remove_redundant_singlar_vertices()
        self.extract_singular_edges()
        # remove some edges
        self.remove_regular_singular_edges()
        self.build_singular_edge_neighbor_relation()

    def is_singular_edge_exist(self, edge_link):
        sid = -1
        for i, es in enumerate(self.singular_es):
            if set(edge_link) == set(es.es_link):
                return i
        return sid

    def add_singular_edges(self, vert_link, edge_link):
        if len(edge_link) == 0:
            return False
        new_seid = self.is_singular_edge_exist(edge_link)
        if new_seid != -1:
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
        return True

    def build_vert_singular_edge_map(self):
        self.vertex_singular_edge_map = {}
        for i, se in enumerate(self.singular_es):
            # build vert and singular edge map
            for v in se.vs_link:
                if v in self.vertex_singular_edge_map:
                    self.vertex_singular_edge_map[v].add(i)
                else:
                    self.vertex_singular_edge_map[v] = {i}

    def add_singular_vertices(self, vids):
        self.singular_vs.update(vids)

    # two end of the edge are singular vertices
    def get_singular_edges(self, eids):
        singular_eids = set()
        for eid in eids:
            edge = self.mesh.Edges[eid]
            if len(check_common_elements(edge.Vids, self.singular_vs)) == 2:
                singular_eids.add(eid)
        return singular_eids

    def is_vert_redundant(self, vid):
        vert = self.mesh.Vertices[vid]
        eids = check_common_elements(vert.neighboring_Eids, self.in_scope_eids)
        singluar_eids = self.get_singular_edges(eids)
        paired_eids = set()
        is_redundant = False
        for seid in singluar_eids:
            if seid in paired_eids:
                continue
            next_edges = self.mesh.get_next_edge(vid, seid)
            next_eid = check_common_elements(self.in_scope_eids, next_edges)
            if len(next_eid) != 1:
                break
            next_eid = next_eid.pop()
            if next_eid in paired_eids:
                continue
            # check valence
            v_0 = check_common_elements(self.mesh.Edges[seid].neighboring_Fids, self.face_scope)
            v_0 = len(v_0)
            v_1 = check_common_elements(self.mesh.Edges[next_eid].neighboring_Fids, self.face_scope)
            v_1 = len(v_1)
            if v_0 == v_1:
                paired_eids.update([seid, next_eid])
        if len(paired_eids) == len(singluar_eids):
            is_redundant = True
        return is_redundant

    def is_vert_at_boundary_corner(self, vid):
        if vid not in self.singular_vs:
            return False
        vert = self.mesh.Vertices[vid]
        nfids = check_common_elements(vert.neighboring_Fids, self.face_scope)
        if nfids == 1:
            return True
        return False

    # edge -> next edge -> remove the vert
    def remove_redundant_singlar_vertices(self):
        print("removing 2D surface redundant singular vertices ...")
        self.redundant_singular_vs = set()
        for vid in self.singular_vs:
           if self.is_vert_redundant(vid):
               self.redundant_singular_vs.add(vid)
           elif self.is_vert_at_boundary_corner(vid):
               self.redundant_singular_vs.add(vid)
        for vid in self.redundant_singular_vs:
            self.singular_vs.remove(vid)
        print("done ...")

    def extract_singular_vertices(self):
        print("extracting 2D surface singular vertices ...")
        # find high valence vertices
        svids = set()
        for vid in self.in_scope_vids:
            is_singular = False
            vert = self.mesh.Vertices[vid]
            fids = check_common_elements(vert.neighboring_Fids, self.face_scope)
            if len(fids) != 4 and (vid not in self.boundary_vids):
                is_singular = True
            elif len(fids) != 2 and (vid in self.boundary_vids):
                is_singular = True

            if not is_singular:
                continue

            # is singularity vert
            svids.add(vid)

        self.add_singular_vertices(svids)
        print("done ...")

    # complete singular vert when traced by a singular edge link
    def extract_singular_edges(self):
        print("extracting 2D surface singular edges ...")
        vid_stack = []
        vid_stack.extend(self.singular_vs)

        self.visited_eids = set()
        self.visited_vids = set()

        while vid_stack:
            vid = vid_stack.pop()
            if vid in self.visited_vids:
                continue
            self.visited_vids.add(vid)
            # if a vert is traced by a edge
            # and the vert is not in singular vs list
            # add it in
            if vid not in self.singular_vs:
                self.redundant_singular_vs.remove(vid)
                self.singular_vs.add(vid)
            vert = self.mesh.Vertices[vid]
            # tracing all edges neighboring with singular vert
            for eid in vert.neighboring_Eids:
                if eid in self.visited_eids:
                    continue
                if eid not in self.in_scope_eids:
                    continue
                vert_link, edge_link = self.trace_edge(vid, eid)
                if not self.add_singular_edges(vert_link, edge_link):
                    continue
                for i in vert_link:
                    if i in self.redundant_singular_vs:
                        vid_stack.append(i)
        self.build_vert_singular_edge_map()
        print("done ...")

    def trace_next_edge(self, vid, eid):
        if (vid in self.singular_vs) or (vid in self.redundant_singular_vs):
            return []
        # next element
        next_edges = self.mesh.get_next_edge(vid, eid)
        next_eid = check_common_elements(self.in_scope_eids, next_edges)
        if len(next_eid) > 1:
            return []
        return next_eid

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

# ----------------------------- remove singular edge
    def get_neighbor_regular_edges(self, vid):
        vert = self.mesh.Vertices[vid]
        # regular mesh eids
        regular_eids = []
        for eid in vert.neighboring_Eids:
            if self.is_edge_regular(eid):
                regular_eids.append(eid)
        return regular_eids

    # only neighboring with 2 quad faces
    def is_edge_regular(self, eid):
        edge = self.mesh.Edges[eid]
        nfids = check_common_elements(edge.neighboring_Fids, self.face_scope)
        if len(nfids) != 2:
            return False
        is_regular = True
        for nfid in nfids:
            face = self.mesh.Faces[nfid]
            if face.is_quad():
                continue
            is_regular = False
            break
        return is_regular

    # remove edges neighboring with non quad face
    def get_next_regular_eids(self, vid, eid):
        if not self.is_edge_regular(eid):
            return []
        next_eids = self.mesh.get_next_edge(vid, eid)
        next_eids = check_common_elements(next_eids, self.in_scope_eids)
        regular_eids = set()
        for eid in next_eids:
            if self.is_edge_regular(eid):
                regular_eids.add(eid)
        return regular_eids

    def get_removable_regular_edges(self, vid):
        neighbor_regular_eids = self.get_neighbor_regular_edges(vid)
        removeable_edge_pairs = []
        for neid in neighbor_regular_eids:
            next_eids = self.get_next_regular_eids(vid, neid)
            if len(next_eids) != 1:
                continue
            edge_pair = {neid}
            edge_pair.update(next_eids)
            if edge_pair in removeable_edge_pairs:
                continue
            removeable_edge_pairs.append(edge_pair)
        return removeable_edge_pairs

    # a singular edge is called regular when the edge
    # contain both non-hex edge and regular edge
    # and after remove non-hex edge, only one regular edge left
    def remove_regular_singular_edges(self):
        print("removing 2D surface redundant singular edges ...")
        # detect candidate vert
        self.vert_and_removeable_regular_edge_pairs = {}
        for vid in self.singular_vs:
            edge_pairs = self.get_removable_regular_edges(vid)
            if len(edge_pairs) > 0:
                self.vert_and_removeable_regular_edge_pairs[vid] = edge_pairs
        # vert - edge pairs
        # to singualr edge pairs
        self.vert_and_removeable_regular_singular_edge_pairs = {}
        for vid, eps in self.vert_and_removeable_regular_edge_pairs.items():
            singular_eids = self.vertex_singular_edge_map[vid]
            new_eps = []
            for ep in eps:
                sep = set()
                for seid in singular_eids:
                    se = self.singular_es[seid]
                    if check_common_elements(se.es_link, ep):
                        sep.add(seid)
                new_eps.append(sep)
            self.vert_and_removeable_regular_singular_edge_pairs[vid] = new_eps

        # check edge chain
        # keep complete singular edge chain when a edge is not removeable
        removeable_singuar_eids = self.check_removeable_regular_singular_edges()
        new_singular_edges = []
        for i, se in enumerate(self.singular_es):
            if i in removeable_singuar_eids:
                continue
            new_singular_edges.append(se)
        self.singular_es = new_singular_edges
        # update vid map
        self.build_vert_singular_edge_map()
        print("done ...")
        
    #
    # 2 is remoables
    # 1 is not
    # 0 is not in the list
    def get_singular_edge_removeable_statue(self, singular_eid):
        se = self.singular_es[singular_eid]
        in_list_vids = check_common_elements(se.start_end_vids,
                                             self.vert_and_removeable_regular_edge_pairs.keys())
        singluar_vids = check_common_elements(se.start_end_vids, self.singular_vs)

        if len(in_list_vids) == 0:
            return 0
        if len(in_list_vids) == 2:
            return 2
        if len(in_list_vids) == 1 and len(singluar_vids) == 1:
            return 2
        return 1

    def check_removeable_regular_singular_edges(self):
        removeable_singular_edges = set()
        unremoveable_singular_edges = set()
        for vid, seps in self.vert_and_removeable_regular_singular_edge_pairs.items():
            for sep in seps:
                for seid in sep:
                    statue = self.get_singular_edge_removeable_statue(seid)
                    if statue == 2:
                        removeable_singular_edges.add(seid)
                        continue
                    if statue == 1:
                        unremoveable_singular_edges.add(seid)

        # find singular edge link
        singular_eid_stack = []
        singular_eid_stack.extend(unremoveable_singular_edges)
        visited_eids = set()
        while singular_eid_stack:
            cur_eids = singular_eid_stack.pop()
            if cur_eids in visited_eids:
                continue
            visited_eids.add(cur_eids)
            edge = self.singular_es[cur_eids]
            for vid in edge.start_end_vids:
                if vid not in self.vert_and_removeable_regular_singular_edge_pairs:
                    continue
                # edge pairs contain the singular edge
                singular_edge_pairs = self.vert_and_removeable_regular_singular_edge_pairs[vid]
                for sep in singular_edge_pairs:
                    # if cur singular edge in the pair
                    # all edges in the pair are not removeable
                    if cur_eids in sep:
                        unremoveable_singular_edges.update(sep)
                        singular_eid_stack.extend(sep)
        # clean eids
        removeable_singular_edges = removeable_singular_edges.difference(unremoveable_singular_edges)
        return removeable_singular_edges

    # ----------------------------- remove singular edge done

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

    def get_valid_neighbor_fids(self, fid):
        face = self.mesh.Faces[fid]
        # valid neighboring fids
        next_fids = set()
        for nfid in face.neighboring_Fids:
            if nfid in self.face_scope:
                continue
            nface = self.mesh.Faces[nfid]
            if check_common_elements(nface.Eids, self.avoid_edge_set):
                continue
            if not check_common_elements(nface.neighboring_Cids, self.cell_scope):
                continue
            next_fids.add(nfid)
        # end
        return next_fids

    # extend face scope to it one ring neighbor
    def extend_one_ring_faces(self):
        for fid in self.input_face_scope:
            next_fids = self.get_valid_neighbor_fids(fid)
            self.extend_faces.update(next_fids)
            self.face_scope.update(next_fids)


    def update_cell_scope(self, cids):
        if not cids:
            return
        self.input_cell_scope = cids
        self.cell_scope = deepcopy(cids)
        self.extend_one_ring_cells()

    def extend_one_ring_cells(self):
        self.extend_cells = set()
        for eid in self.avoid_edge_set:
            edge = self.mesh.Edges[eid]
            next_cids = set(edge.neighboring_Cids).difference(self.cell_scope)
            self.extend_cells.update(next_cids)
            self.cell_scope.update(next_cids)

