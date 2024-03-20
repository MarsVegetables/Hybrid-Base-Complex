from structure_analysis.region_cell_tracer import RegionCellTracer
from mesh_structure.mesh_common import check_common_elements, reduce_id_groups
from mesh_structure.mesh_elements import Interior_Regular_E, Boundary_Regular_E
from base_complex.base_complex_elements import SingularityEdge
import random
from structure_analysis.singular_parallel_edge_grouping import SingularParallelEdgeGrouping


class ValenceSingularityGraph3D:
    def __init__(self, mesh, parallel_eids=None, region_cids=None):
        self.mesh = mesh
        if parallel_eids:
            self.parallel_eids = parallel_eids
        else:
            self.parallel_eids = set()
        # fids get from a shee
        # self.sheet_fids = sheet_fids
        self.region_cids = set()
        self.extended_cids = set()
        if region_cids:
            self.region_cids = set(region_cids)
            if self.parallel_eids:
                extend_region = RegionCellTracer(self.mesh, self.parallel_eids, self.region_cids)
                self.extended_cids = extend_region.extend_non_hex_cells_one_ring()
            # add extended cell into region cell set
            self.region_cids.update(self.extended_cids)

        # region elements
        self.region_fids = set()
        self.region_eids = set()
        self.region_vids = set()
        self.find_region_elements()

        # region boundary elements
        self.boundary_fids = set()
        self.boundary_eids = set()
        self.boundary_vids = set()
        self.find_boundary_elements()

        #
        self.visited_eids = set()
        self.visited_fids = set()
        self.patched_fids = set()

        # singularity graph elements
        # partial singularity edge link are seed singular edges
        self.seed_singular_edges = set()
        # all mesh singular edges
        self.singular_mesh_edges = set()
        self.singular_es = []
        self.singular_vs = set()

        # visualization
        # self.color_num = 0
        # edge color number
        self.edge_color_num = 0
        self.vertex_singular_edge_map = {}
        # mesh edge id to singular edge id
        self.edge_singular_edge_map = {}
        self.vert_and_removable_singular_edge_pairs = {}
        self.vert_and_removable_edge_pairs = {}

        # clean graph
        self.redundant_singular_vs = set()
        self.redundant_singular_es = []
        # the singular edges only contain parallel edges
        self.sheet_parallel_singular_es = []

        # final update flag
        self.final_update_flag = False
        self.singular_edge_group_id = {}

    def extract_valence_singular_graph(self, clean_graph=False, remove_parallel_edges=True):
        self.init_structure_singularity_graph()
        self.complete_valence_singularity_graph()
        if clean_graph:
            self.clean_valence_singularity_graph(remove_parallel_edges)
            self.final_singular_elements()

    def init_structure_singularity_graph(self):
        self.find_seed_singular_edges()
        self.update_singular_edge()
        self.init_singular_vertices()

    def clean_valence_singularity_graph(self, remove_parallel_edges=True):
        self.detect_redundant_singular_vert()
        self.remove_redundant_singular_edges()
        if remove_parallel_edges:
            self.remove_sheet_parallel_singular_edges()
        self.update_mesh_singular_edge_set()
        # update vid map
        self.init_singular_vertices()
        self.build_vert_singular_edge_map()
        self.build_singular_edge_neighbor_relation()

    def clean_singular_edge_neighbor_info(self):
        for se in self.singular_es:
            se.init_neighbors()

    def update_singular_edge(self):
        self.singular_es = []
        self.extract_singular_edges(self.seed_singular_edges)
        self.build_vert_singular_edge_map()
        self.build_singular_edge_neighbor_relation()

    def final_singular_elements(self):
        self.final_update_flag = True
        self.update_final_singular_v()
        self.update_seed_singular_edges()
        self.update_singular_edge()
        self.build_vert_singular_edge_map()
        self.build_singular_edge_neighbor_relation()

    def find_region_elements(self):
        # region is entire mesh
        if not self.region_cids:
            self.region_cids = self.mesh.Cells.keys()
            self.region_fids = self.mesh.Faces.keys()
            self.region_eids = self.mesh.Edges.keys()
            self.region_vids = self.mesh.Vertices.keys()
            return
        for cid in self.region_cids:
            cell = self.mesh.Cells[cid]
            self.region_vids.update(cell.Vids)
            self.region_eids.update(cell.Eids)
            self.region_fids.update(cell.Fids)

    def find_boundary_elements(self):
        for fid in self.region_fids:
            face = self.mesh.Faces[fid]
            ncids = self.get_neighboring_cids(face)
            if len(ncids) != 2:
                self.boundary_fids.add(fid)
                self.boundary_vids.update(face.Vids)
                self.boundary_eids.update(face.Eids)

    # def complete_valence_singularity_graph(self):
    #     print("completing valence singularity graph ...")
    #     all_traced = False
    #     while not all_traced:
    #         # trace all edge sround singular vertex
    #         self.complete_seed_edges()
    #         # update new singular vertices
    #         self.complete_singular_vert()
    #         all_traced = self.check_edge_surround_singular_vertices()
    #     self.update_singular_edge()
    #     print("done ...")

# ------------------------------------------------------------------------------------------------------------------ complete valence SG
    def complete_valence_singularity_graph(self):
        print("completing valence singularity graph ...")
        # trace seed edges by trace faces
        self.build_face_patches()
        # find all seed edges
        self.build_patch_eids()
        # update new singular vertices
        self.complete_singular_vert()
        self.complete_seed_edges()
        self.update_singular_edge()
        print("done ...")

    def build_face_patches(self):
        print("extracting VSG face patches ...")
        self.visited_eids = set()
        self.visited_fids = set()
        self.face_patch_fids = [] # for debug
        self.separated_face_patches = []
        start_eids = set(self.singular_mesh_edges)
        for vid in self.singular_vs:
            vert = self.mesh.Vertices[vid]
            neids = self.get_neighboring_eids(vert)
            start_eids.update(neids)
        for eid in start_eids:
            self.visited_eids.add(eid)
            edge = self.mesh.Edges[eid]
            nfids = check_common_elements(edge.neighboring_Fids, self.region_fids)
            for fid in nfids:
                if fid in self.visited_fids:
                    continue
                pathFids = self.trace_patch_faces(fid)
                self.face_patch_fids.append(pathFids) # for debug

        for fid in self.region_fids:
            if (fid in self.visited_fids) or (fid in self.boundary_fids):
                self.patched_fids.add(fid)
        print("done ...")

    def trace_patch_faces(self, startFid):
        fidStack = [startFid]
        patchFids = set(fidStack)
        while len(fidStack) > 0:
            cur_fid = fidStack.pop()
            if cur_fid in self.visited_fids:
                continue
            self.visited_fids.add(cur_fid)
            op_fids = self.mesh.get_opposite_faces_of_face(cur_fid)
            op_fids = check_common_elements(op_fids, self.region_fids)
            if len(op_fids) > 0:
                fidStack.extend(op_fids)
                patchFids.update(op_fids)
        return patchFids

    def get_opposite_faces_of_face(self, fid):
        oppositeFaces = []
        face = self.Faces[fid]
        for eid in face.Eids:
            oppositeFid = self.get_opposite_face_on_regular_eid(fid, eid)
            if len(oppositeFid) > 0:
                oppositeFaces.extend(oppositeFid)
        return oppositeFaces

    def get_opposite_face_on_regular_eid(self, fid, eid):
        # singularity edge should only contain one opposite face of an edge at a time
        oppositeFace = []
        edge = self.Edges[eid]
        if (eid in self.singular_mesh_edges) or (eid in self.boundary_eids):
            return oppositeFace
        oppositeFace = self.get_opposite_face_on_any_type_eid(fid, eid)
        # only contain one opposite face
        return oppositeFace

    def get_opposite_face_on_any_type_eid(self, fid, eid):
        # singularity edge should only contain one opposite face of an edge at a time
        oppositeFace = []
        face = self.Faces[fid]
        edge = self.Edges[eid]
        nfids = check_common_elements(edge.neighboring_Fids, self.region_fids)
        for nfid in nfids:
            if nfid in self.visited_fids:
                continue
            nface = self.Faces[nfid]
            sharingCids = check_common_elements(face.neighboring_Cids, nface.neighboring_Cids)
            sharingCids = check_common_elements(sharingCids, self.region_cids)
            # if both face and nface contain same neighboring cids
            # they are at the same cid
            if len(sharingCids) > 0:
                continue
            oppositeFace.append(nfid)
            # if I am sure that a regular edge only has one opposite face
            # I can return here
            # return nfid
        if len(oppositeFace) > 1:
            pass
            # print("a regular edge contain a face have two opposite face....")
        # only contain one opposite face
        return oppositeFace

    # https://github.com/Cotrik/CotrikMesh/blob/d890c9e596fa3242342b8a1b7fc7b9ce194f03a7/libcotrik/src/BaseComplex.cpp#L1420
    def build_patch_eids(self):
        print("extracting complete seed edges ...")
        if not self.patched_fids:
            print("No Fids is added into patch ...")
            return
        self.visited_eids = set()
        tmp_eids = set()
        #
        for fid in self.patched_fids:
            face = self.mesh.Faces[fid]
            for eid in face.Eids:
                tmp_eids.add(eid)
        remove_eids = set()
        for eid in tmp_eids:
            edge = self.mesh.Edges[eid]
            nfids = check_common_elements(edge.neighboring_Fids, self.region_fids)
            fid_dif = nfids - self.patched_fids
            if fid_dif:
                remove_eids.add(eid)
        # remove
        for e in remove_eids:
            tmp_eids.remove(e)
        self.seed_singular_edges = list(tmp_eids)
        print("done ...")

    # # trace other edges from singular vertices
    def complete_seed_edges(self):
        self.seed_singular_edges = set()
        self.visited_eids = set()
        # edge surround singular vertices
        for vid in self.singular_vs:
            vert = self.mesh.Vertices[vid]
            vert_neids = self.get_neighboring_eids(vert)
            for eid in vert_neids:
                if eid not in self.region_eids:
                    continue
                self.seed_singular_edges.add(eid)

# ----------------------------------------------------------------------------------------------- complete VSG done

    # ------------------------------ singularity graph related
    def check_edge_surround_singular_vertices(self):
        all_traced = True
        for vid in self.singular_vs:
            if vid not in self.region_vids:
                continue
            vert = self.mesh.Vertices[vid]
            in_scope_eids = check_common_elements(self.region_eids, vert.neighboring_Eids)
            tmp = check_common_elements(in_scope_eids, self.seed_singular_edges)
            if in_scope_eids == tmp:
                continue
            all_traced = False
        return all_traced

    # only consider edge valence
    # edge neighboring cell type is not considered
    def find_seed_singular_edges(self):
        self.seed_singular_edges = set()
        for eid in self.region_eids:
            edge = self.mesh.Edges[eid]
            valence = len(check_common_elements(edge.neighboring_Cids, self.region_cids))
            if (eid in self.boundary_eids) and (valence != Boundary_Regular_E):
                self.seed_singular_edges.add(eid)
            elif (eid not in self.boundary_eids) and (valence != Interior_Regular_E):
                self.seed_singular_edges.add(eid)

    # in 3d space, if the edge contain more than one next edge
    # the edge is singularity
    def trace_next_edge(self, vid, eid):
        if eid not in self.region_eids:
            return []
        if vid not in self.region_vids:
            return []
        if vid in self.singular_vs:
            return []
        # next element
        next_eid = self.mesh.get_next_edge(vid, eid)
        next_eid = check_common_elements(next_eid, self.region_eids)
        if self.final_update_flag:
            new_next_eids = set()
            for neid in next_eid:
                if self.valid_final_mesh_edge(neid):
                    new_next_eids.add(neid)
            next_eid = new_next_eids
        if len(next_eid) != 1:
            return []
        if self.final_update_flag:
            return next_eid
        # same valence
        edge = self.mesh.Edges[eid]
        scope_cid_1 = check_common_elements(edge.neighboring_Cids, self.region_cids)
        next_eid = next_eid.pop()
        nedge = self.mesh.Edges[next_eid]
        scope_cid_2 = check_common_elements(nedge.neighboring_Cids, self.region_cids)
        if len(scope_cid_1) != len(scope_cid_2):
            return []
        return [next_eid]

    # give an edge, trace to both direction
    def trace_edge(self, vid, eid):
        # sequence of edge is matter
        edge_link = []
        vert_link = [vid]
        next_eid = [eid]
        next_vid = vid
        while next_eid:
            cur_eid = next_eid.pop()
            if cur_eid not in self.region_eids:
                continue
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
            next_eid.extend(self.trace_next_edge(next_vid, cur_eid))
        return vert_link, edge_link

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

    def extract_singular_edges(self, eids):
        print("extracting singular edges ...")
        self.visited_eids = set()
        for eid in eids:
            if eid not in self.region_eids:
                continue
            if eid in self.visited_eids:
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
        self.update_mesh_singular_edge_set()
        print("done ...")

    def init_singular_vertices(self):
        self.singular_vs = set()
        for se in self.singular_es:
            self.singular_vs.update(se.start_end_vids)
        # for v in self.singular_vs:
        #     if v in self.redundant_singular_vs:
        #         self.redundant_singular_vs.remove(v)

    def is_two_edge_opposite(self, e0, e1):
        edge_0 = self.mesh.Edges[e0]
        edge_1 = self.mesh.Edges[e1]
        common_fids = check_common_elements(edge_0.neighboring_Fids, edge_1.neighboring_Fids)
        common_fids = check_common_elements(common_fids, self.region_fids)
        if common_fids:
            return False
        return True

    def complete_singular_vert(self):
        self.singular_vs = set()
        for eid in self.seed_singular_edges:
            edge = self.mesh.Edges[eid]
            for vid in edge.Vids:
                vert = self.mesh.Vertices[vid]
                in_scope_eids = check_common_elements(vert.neighboring_Eids, self.region_eids)
                nse = check_common_elements(in_scope_eids, self.seed_singular_edges)
                if (len(nse) > 0) and (len(nse) != 2):
                    self.singular_vs.add(vid)
                    continue
                if (len(nse) == 2) and (not self.is_two_edge_opposite(nse.pop(), nse.pop())):
                    self.singular_vs.add(vid)
                    continue


    def build_vert_singular_edge_map(self):
        self.vertex_singular_edge_map = {}
        for i, se in enumerate(self.singular_es):
            # build vert and singular edge map
            for v in se.vs_link:
                if v in self.vertex_singular_edge_map:
                    self.vertex_singular_edge_map[v].add(i)
                else:
                    self.vertex_singular_edge_map[v] = {i}

    # assign edge color
    def build_singular_edge_neighbor_relation(self):
        self.clean_singular_edge_neighbor_info()
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
        self.edge_color_num = color_number
        return color_number

    def assign_singular_edge_color(self):
        print("assigning singularity edge color ...")
        if self.parallel_eids:
            self.assign_singular_edge_color_group()
        else:
            self.assign_singular_edge_color_random()
        print("assigning singularity edge color ... done")

    def assign_singular_edge_color_random(self):
        for se in self.singular_es:
            neighbor_colors = set()
            for n_se_id in se.neighboring_Eids:
                cc = self.singular_es[n_se_id]
                neighbor_colors.add(cc.color)
            candidate_colors = set()
            for i in range(self.edge_color_num):
                # keep rad color reserved
                if i == 0:
                    continue
                if i not in neighbor_colors:
                    candidate_colors.add(i)
            picked_color = random.choice(list(candidate_colors))
            se.color = picked_color

    def assign_singular_edge_color_group(self):
        singular_edge_grouper = SingularParallelEdgeGrouping(self)
        singular_edge_grouper.grouping(0.8)
        self.singular_edge_group_id = singular_edge_grouper.singular_edge_group_id
        group_color = {}
        for i, se in enumerate(self.singular_es):
            neighbor_colors = set()
            for n_se_id in se.neighboring_Eids:
                cc = self.singular_es[n_se_id]
                neighbor_colors.add(cc.color)
            candidate_colors = set()

            # group color
            new_group_color = False
            gid = self.singular_edge_group_id[i]
            if gid in group_color:
                candidate_colors.add(group_color[gid])

            # color
            if not candidate_colors:
                new_group_color = True
                for j in range(self.edge_color_num):
                    # keep rad color reserved
                    if j == 0:
                        continue
                    if j not in neighbor_colors:
                        candidate_colors.add(j)
            picked_color = random.choice(list(candidate_colors))
            se.color = picked_color

            # add color to group
            if new_group_color:
                gid = self.singular_edge_group_id[i]
                group_color[gid] = picked_color

    def is_vert_at_corner(self, vid):
        vert = self.mesh.Vertices[vid]
        vert_ncids = self.get_neighboring_cids(vert)
        if len(vert_ncids) == 1:
            return True
        return False

    # clean graph
    def is_vert_at_hex_corner(self, vid):
        if vid not in self.singular_vs:
            return False
        vert = self.mesh.Vertices[vid]
        ncids = check_common_elements(vert.neighboring_Cids, self.region_cids)
        if len(ncids) == 1:
            ncid = ncids.pop()
            if self.mesh.Cells[ncid].is_hexa():
                return True
        return False

    def is_vert_at_all_quad_corner(self, vid):
        if vid not in self.singular_vs:
            return False
        vert = self.mesh.Vertices[vid]
        ncids = check_common_elements(vert.neighboring_Cids, self.region_cids)
        if len(ncids) != 1:
            return False
        nq = 0
        nfids = check_common_elements(vert.neighboring_Fids, self.region_fids)
        for nf in nfids:
            nface = self.mesh.Faces[nf]
            if nface.is_quad():
                nq += 1
        nqr = nq / len(nfids)
        if nqr != 1:
            return False
        return True

    # detect edge neighbor cell types
    def is_singular_edge_regular_valence(self, eid):
        if eid not in self.singular_mesh_edges:
            return False
        edge = self.mesh.Edges[eid]
        ncids = check_common_elements(edge.neighboring_Cids, self.region_cids)
        if eid in self.boundary_eids and len(ncids) == Boundary_Regular_E:
            return True
        if eid not in self.boundary_eids and len(ncids) == Interior_Regular_E:
            return True
        return False

    def is_edge_in_bnd_corner_with_all_quad_neighbor(self, eid):
        edge = self.mesh.Edges[eid]
        ncids = check_common_elements(edge.neighboring_Cids, self.region_cids)
        if len(ncids) != 1:
            return False
        nqr = self.get_neighbor_quad_ratio(eid)
        if nqr != 1:
            return False
        return True

    def vert_only_neighbor_with_regular_valance_singular_edges(self, vid):
        inregular_valence_edge = set()
        vert = self.mesh.Vertices[vid]
        mesh_singular_edges = check_common_elements(vert.neighboring_Eids, self.singular_mesh_edges)
        for eid in mesh_singular_edges:
            if self.is_singular_edge_regular_valence(eid):
                continue
            if self.is_edge_in_bnd_corner_with_all_quad_neighbor(eid):
                continue
            inregular_valence_edge.add(eid)
        if len(inregular_valence_edge) == 0:
            return True
        return False

    def detect_redundant_singular_vert(self):
        print("detecting redundant singular vertex ...")
        self.redundant_singular_vs = set()
        for vid in self.singular_vs:
            # if (self.is_vert_at_hex_corner(vid) or
            #         self.vert_only_neighbor_with_regular_valance_singular_edges(vid)):
            if (self.is_vert_at_all_quad_corner(vid) or
                    self.vert_only_neighbor_with_regular_valance_singular_edges(vid)):
                self.redundant_singular_vs.add(vid)
        for vid in self.redundant_singular_vs:
            self.singular_vs.remove(vid)
        print("done ...")

    def update_mesh_singular_edge_set(self):
        self.singular_mesh_edges = set()
        self.edge_singular_edge_map = {}
        for i, se in enumerate(self.singular_es):
            self.singular_mesh_edges.update(se.es_link)
            for eid in se.es_link:
                if eid not in self.edge_singular_edge_map:
                    self.edge_singular_edge_map[eid] = set()
                self.edge_singular_edge_map[eid].add(i)

    # ----------------- redundant singular edge

    def get_neighbor_singular_edges_with_regular_valence(self, vid):
        vert = self.mesh.Vertices[vid]
        # regular mesh eids
        regular_eids = []
        vert_neid = self.get_neighboring_eids(vert)
        for eid in vert_neid:
            if self.is_singular_edge_regular_valence(eid):
                regular_eids.append(eid)
            # or in hex cell corner
            if self.is_edge_in_bnd_corner_with_all_quad_neighbor(eid):
                regular_eids.append(eid)
        return regular_eids

    # remove edges neighboring with non quad face
    def get_next_region_eids(self, vid, eid):
        # if not self.is_edge_regular_and_all_hex_neighbor(eid):
        #     return []
        next_eids = self.mesh.get_next_edge(vid, eid)
        next_eids = check_common_elements(next_eids, self.region_eids)
        return next_eids

    def is_edge_all_hex_neighbor(self, eid):
        edge = self.mesh.Edges[eid]
        all_hex = True
        edge_ncids = self.get_neighboring_cids(edge)
        for ncid in edge_ncids:
            cell = self.mesh.Cells[ncid]
            if cell.is_hexa():
                continue
            all_hex = False
        return all_hex

    def is_edge_all_quad_neighbor(self, eid):
        nqr = self.get_neighbor_quad_ratio(eid)
        if nqr == 1:
            return True
        return False

    # get_neighbor_regular_valence_singular_edge_pairs
    # def get_candidate_removable_edge_pairs(self, vid):
    #     neighbor_regular_eids = self.get_neighbor_singular_edges_with_regular_valence(vid)
    #     removable_edge_pairs = []
    #     for neid in neighbor_regular_eids:
    #         if self.get_neighbor_quad_ratio(neid) != 1:
    #             continue
    #         next_eids = self.get_next_region_eids(vid, neid)
    #         next_hex_eids = set()
    #         # if multiple next edges
    #         if len(next_eids) != 1:
    #             # find next all hex neighbor edge
    #             for e in next_eids:
    #                 # if self.is_edge_all_hex_neighbor(e):
    #                 if self.is_edge_all_quad_neighbor(e):
    #                     next_hex_eids.add(e)
    #         if len(next_eids) != 1 and len(next_hex_eids) != 1:
    #             continue
    #         # if only one next all hex neighbor edge
    #         if len(next_eids) != 1:
    #             next_eids = next_hex_eids
    #         # the edge valence must regular
    #         if not check_common_elements(next_eids, neighbor_regular_eids):
    #             continue
    #         edge_pair = {neid}
    #         edge_pair.update(next_eids)
    #         if edge_pair in removable_edge_pairs:
    #             continue
    #         removable_edge_pairs.append(edge_pair)
    #     return removable_edge_pairs

    # all next edge must be all quad neighbor
    def get_candidate_removable_edge_pairs(self, vid):
        neighbor_regular_eids = self.get_neighbor_singular_edges_with_regular_valence(vid)
        removable_edge_pairs = []
        for neid in neighbor_regular_eids:
            if self.get_neighbor_quad_ratio(neid) != 1:
                continue
            next_eids = self.get_next_region_eids(vid, neid)
            next_hex_eids = set()
            # if multiple next edges
            if len(next_eids) != 0:
                # find next all hex neighbor edge
                for e in next_eids:
                    # if self.is_edge_all_hex_neighbor(e):
                    if self.is_edge_all_quad_neighbor(e):
                        next_hex_eids.add(e)
            if len(next_eids) != 1 and len(next_hex_eids) != 1:
                continue
            # if only one next all hex neighbor edge
            if len(next_eids) != 1:
                next_eids = next_hex_eids
            # the edge valence must regular
            if not check_common_elements(next_eids, neighbor_regular_eids):
                continue
            edge_pair = {neid}
            edge_pair.update(next_eids)
            if edge_pair in removable_edge_pairs:
                continue
            removable_edge_pairs.append(edge_pair)
        return removable_edge_pairs

    def get_neighbor_quad_ratio(self, eid):
        edge = self.mesh.Edges[eid]
        nq = 0
        edge_nfids = self.get_neighboring_fids(edge)
        for f in edge_nfids:
            face = self.mesh.Faces[f]
            if face.is_quad():
                nq += 1
        return nq / len(edge_nfids)

    def get_neighbor_quad_num(self, eid):
        edge = self.mesh.Edges[eid]
        nq = 0
        edge_nfids = self.get_neighboring_fids(edge)
        for f in edge_nfids:
            face = self.mesh.Faces[f]
            if face.is_quad():
                nq += 1
        return nq

    # edge
    # gap is 0
    # if valence 4 internal
    # if valence 2 bnd
    # if valnce 1 corner
    def get_edge_valence_gap(self, eid):
        gap = 0
        edge = self.mesh.Edges[eid]
        ncids = len(check_common_elements(edge.neighboring_Cids, self.region_cids))
        if edge.is_boundary:
            if ncids == 1:
                return gap
            return abs(ncids - Boundary_Regular_E)
        return abs(ncids - Interior_Regular_E)

    def get_edge_pair_quad_ratio_and_valence_gap(self, edge_paris):
        quad_ratio = 0
        valence_gap = 0
        for e in edge_paris:
            valence_gap += self.get_edge_valence_gap(e)
            nqr = self.get_neighbor_quad_ratio(e)
            quad_ratio += nqr
        return quad_ratio, valence_gap

    def verify_candidate_removable_edge_pairs(self):
        new_v_e_pairs = {}
        for v, eps in self.vert_and_removable_edge_pairs.items():
            grouped_eps = reduce_id_groups(eps)
            new_eps = []
            for gep in grouped_eps:
                if len(gep) == 2:
                    new_eps.append(gep)
                    continue
                if len(gep) < 2:
                    continue
                in_group_pairs = []
                for ep in eps:
                    if check_common_elements(ep, gep):
                        in_group_pairs.append(ep)
                # pair valence gap, and quad ratio
                quad_ratio_list = []
                valence_gap_list = []
                for ep in in_group_pairs:
                    qr, vg = self.get_edge_pair_quad_ratio_and_valence_gap(ep)
                    quad_ratio_list.append(qr)
                    valence_gap_list.append(vg)
                # lowest valence gap
                pids = set()
                sort_vg = sorted(valence_gap_list)
                for i, vg in enumerate(valence_gap_list):
                    if vg == sort_vg[0]:
                        pids.add(i)
                if len(pids) == 1:
                    new_eps.append(in_group_pairs[pids.pop()])
                    continue

                # update quad ratio list
                pair_ids_list = []
                remain_quad_ratio_list = []
                for i, qr in enumerate(quad_ratio_list):
                    if i not in pids:
                        continue
                    pair_ids_list.append(i)
                    remain_quad_ratio_list.append(qr)
                # highest quad ratio pair id
                sort_qr = sorted(remain_quad_ratio_list, reverse=True)
                pids = set()
                for i, qr in enumerate(remain_quad_ratio_list):
                    if qr == sort_qr[0]:
                        pid = pair_ids_list[i]
                        pids.add(pid)
                if len(pids) != 1:
                    continue
                pid = pids.pop()
                new_eps.append(in_group_pairs[pid])
            if len(new_eps) == 0:
                continue
            new_v_e_pairs[v] = new_eps
        self.vert_and_removable_edge_pairs = new_v_e_pairs

    def remove_redundant_singular_edges(self):
        print("removing 3D redundant singular edges ...")
        # detect candidate vert
        # mesh level singular edges
        # each eid in the dict are mesh edge id
        # vert_and_removable_regular_valence_singular_edge_pairs
        self.vert_and_removable_edge_pairs = {}
        for vid in self.singular_vs:
            edge_pairs = self.get_candidate_removable_edge_pairs(vid)
            if len(edge_pairs) > 0:
                self.vert_and_removable_edge_pairs[vid] = edge_pairs

        # check the candidate pairs
        self.verify_candidate_removable_edge_pairs()

        # vert - edge pairs
        # to singular edge pairs
        # se id in self.singular_es
        self.vert_and_removable_singular_edge_pairs = {}
        for vid, eps in self.vert_and_removable_edge_pairs.items():
            singular_eids = self.vertex_singular_edge_map[vid]
            new_eps = []
            for ep in eps:
                sep = set()
                for seid in singular_eids:
                    se = self.singular_es[seid]
                    if check_common_elements(se.es_link, ep):
                        sep.add(seid)
                new_eps.append(sep)
            self.vert_and_removable_singular_edge_pairs[vid] = new_eps

        # check edge chain
        # keep complete singular edge chain when an edge is not removable
        removable_singular_eids = self.check_removable_singular_edges()
        new_singular_edges = []
        self.redundant_singular_es = []
        for i, se in enumerate(self.singular_es):
            if i in removable_singular_eids:
                self.redundant_singular_es.append(se)
                continue
            new_singular_edges.append(se)
        self.singular_es = new_singular_edges
        print("done ...")

    def get_singular_edge_redandent_vert_num(self, singular_eid):
        se = self.singular_es[singular_eid]
        # edge two ends are redundant vert
        in_list_vids = check_common_elements(se.start_end_vids,
                                             self.redundant_singular_vs)
        return len(in_list_vids)

    def mesh_eids_to_singular_eids(self, mesh_eids):
        singular_eids = set()
        for eid in mesh_eids:
            if eid not in self.edge_singular_edge_map:
                continue
            s = self.edge_singular_edge_map[eid]
            if len(s) != 1:
                print(f"mesh edge {eid} contain {len(s)} singular edges : {s}")
            singular_eids.update(s)
        return singular_eids

    # two edge, neighboring with each other, and sharing cells
    # but does not share any boundary faces
    # works only when the singular edge does not have any valid next singular edge at the input vert
    def get_same_neighbor_cell_bnd_singular_eid(self, vid, singular_eid):
        if vid not in self.boundary_vids:
            return []
        if not self.singular_es[singular_eid].is_boundary:
            return []
        # both vid and singular edge at boundary
        vert = self.mesh.Vertices[vid]
        vert_neids = self.get_neighboring_eids(vert)
        se = self.singular_es[singular_eid]
        input_mesh_edge = check_common_elements(se.es_link, vert_neids)
        if len(input_mesh_edge) != 1:
            print(f"vert {vid} has {len(input_mesh_edge)} neighboring with singular edge {singular_eid}")
        input_mesh_eid = input_mesh_edge.pop()
        input_mesh_edge = self.mesh.Edges[input_mesh_eid]
        input_mesh_edge_ncids = self.get_neighboring_cids(input_mesh_edge)
        next_mesh_eids = set()
        next_parallel_eids = set()
        for eid in vert_neids:
            if eid == input_mesh_eid:
                continue
            if eid not in self.boundary_eids:
                continue
            tmp_edge = self.mesh.Edges[eid]
            tmp_ncids = self.get_neighboring_cids(tmp_edge)
            if set(input_mesh_edge_ncids) != set(tmp_ncids):
                continue
            next_mesh_eids.add(eid)
            if eid in self.parallel_eids:
                next_parallel_eids.add(eid)
        # remove parallel edge from corner trace
        if len(next_mesh_eids) > 1:
            next_mesh_eids = next_mesh_eids.difference(next_parallel_eids)
        # mesh eids to singular edges
        next_seids = self.mesh_eids_to_singular_eids(next_mesh_eids)
        return next_seids

    def get_next_singular_eid_one_side(self, vid, singular_eid):
        mesh_eid = set()
        se = self.singular_es[singular_eid]
        for eid in se.es_link:
            edge = self.mesh.Edges[eid]
            if vid in edge.Vids:
                mesh_eid.add(eid)
        if len(mesh_eid) != 1:
            print(f"se {singular_eid} contain edge {mesh_eid} connected to vert {vid}")
        mesh_eid = mesh_eid.pop()
        # all next mesh edges
        next_eids = self.get_next_region_eids(vid, mesh_eid)
        if not next_eids:
            return set()
        next_singular_eids = self.mesh_eids_to_singular_eids(next_eids)
        #
        if vid not in self.vert_and_removable_singular_edge_pairs:
            return next_singular_eids
        # edge contain paris on the end
        related_edge_pairs = []
        paired_singular_eids = set()
        # find realted singular edge pairs
        pairs = self.vert_and_removable_singular_edge_pairs[vid]
        related_edge_pairs.extend(pairs)
        for ep in pairs:
            paired_singular_eids.update(ep)
        # paired edge only been active by its pair
        # if no paired edge is found
        # active not in pair edges
        not_paired_eids = set()
        for seid in next_singular_eids:
            if seid in paired_singular_eids:
                continue
            not_paired_eids.add(seid)
        valid_paired_eids = set()
        for seid in next_singular_eids:
            pair = {singular_eid, seid}
            if pair in related_edge_pairs:
                valid_paired_eids.add(seid)
        if valid_paired_eids:
            return valid_paired_eids
        return not_paired_eids

    def get_next_singular_eid(self, singular_eid):
        se = self.singular_es[singular_eid]
        next_singular_edges = set()
        # find related singular edge pairs
        for vid in se.start_end_vids:
            next_eids = self.get_next_singular_eid_one_side(vid, singular_eid)
            if len(next_eids) == 0:
                next_eids = self.get_same_neighbor_cell_bnd_singular_eid(vid, singular_eid)
            next_singular_edges.update(next_eids)
        return next_singular_edges

    def check_removable_singular_edges(self):
        removable_singular_edges = set()
        unremovable_singular_edges = set()
        # edge contain valid next candidate removable edge
        for vid, seps in self.vert_and_removable_singular_edge_pairs.items():
            for ep in seps:
                removable_singular_edges.update(ep)
        for i, se in enumerate(self.singular_es):
            if i in removable_singular_edges:
                continue
            rv = self.get_singular_edge_redandent_vert_num(i)
            if rv == 2:
                removable_singular_edges.add(i)
                continue
            unremovable_singular_edges.add(i)

        self.removable_singular_edge_objs = [self.singular_es[i] for i in removable_singular_edges]
        self.unremovable_singular_edges_objs = [self.singular_es[i] for i in unremovable_singular_edges]
        # find singular edge link
        singular_eid_stack = []
        singular_eid_stack.extend(unremovable_singular_edges)
        # new unremovable singular edges
        unremovable_singular_edges = set()
        while singular_eid_stack:
            cur_eid = singular_eid_stack.pop()
            if cur_eid in unremovable_singular_edges:
                continue
            unremovable_singular_edges.add(cur_eid)
            next_singular_eids = self.get_next_singular_eid(cur_eid)
            singular_eid_stack.extend(next_singular_eids)
        # clean eids
        removable_singular_edges = removable_singular_edges.difference(unremovable_singular_edges)
        return removable_singular_edges

    def remove_sheet_parallel_singular_edges(self):
        if not self.parallel_eids:
            return
        new_se = []
        remove_id = set()
        for i, se in enumerate(self.singular_es):
            tmp = check_common_elements(se.es_link, self.parallel_eids)
            if tmp == set(se.es_link):
                self.sheet_parallel_singular_es.append(se)
                remove_id.add(i)

        for i, se in enumerate(self.singular_es):
            if i in remove_id:
                continue
            new_se.append(se)
        self.singular_es = new_se

    def detect_valid_singular_v(self, vid):
        nseid = self.vertex_singular_edge_map[vid]
        if len(nseid) != 2:
            return True
        vert = self.mesh.Vertices[vid]
        used_neids = check_common_elements(self.singular_mesh_edges, vert.neighboring_Eids)
        for ncid in vert.neighboring_Cids:
            cell = self.mesh.Cells[ncid]
            common_eids = check_common_elements(cell.Eids, used_neids)
            if len(common_eids) > 1:
                return True
        return False

    def update_final_singular_v(self):
        # self.build_vert_singular_edge_map()
        new_singular_v = set()
        for v in self.singular_vs:
            # nseid = self.vertex_singular_edge_map[v]
            if self.detect_valid_singular_v(v):
                new_singular_v.add(v)
        self.singular_vs = new_singular_v

    def valid_final_mesh_edge(self, eid):
        if eid not in self.region_eids:
            return False
        if not self.final_update_flag:
            return False
        if eid in self.singular_mesh_edges:
            return True
        return False

    def update_seed_singular_edges(self):
        self.seed_singular_edges = set()
        for v in self.singular_vs:
            vert = self.mesh.Vertices[v]
            vert_neids = self.get_neighboring_eids(vert)
            seeds = check_common_elements(vert_neids, self.singular_mesh_edges)
            self.seed_singular_edges.update(seeds)

    def get_neighboring_cids(self, element):
        ncids = check_common_elements(element.neighboring_Cids, self.region_cids)
        return ncids

    def get_neighboring_fids(self, element):
        nfids = check_common_elements(element.neighboring_Fids, self.region_fids)
        return nfids

    def get_neighboring_eids(self, element):
        neids = check_common_elements(element.neighboring_Eids, self.region_eids)
        return neids

    def get_neighboring_vids(self, element):
        nvids = check_common_elements(element.neighboring_Vids, self.region_vids)
        return nvids



