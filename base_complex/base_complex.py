import random

from mesh_structure.mesh_common import check_common_elements
from base_complex.singularity_structure import SingularityStructure
from base_complex.base_complex_elements import ComponentCell, ComponentFace, ComponentEdge, ComponentVert


class BaseComplex:
    def __init__(self, singularity_structure, iteration_num=0):
        self.singularity_structure = singularity_structure
        self.mesh = self.singularity_structure.mesh
        self.init_attr()
        self.iter_num = iteration_num
        print("===================== iteration : {} =========================".format(self.iter_num))

    def init_attr(self):
        self.Fids = []  # base complex Fids, mesh face id
        self.Eids = []  # base complex Eids, mesh edge id
        self.Vids = []  # base complex Vids, mesh vertex id
        self.Cids = []  # base complex Cids, mesh cell id
        self.componentCs = []  # base complex component cells, component cell id
        self.componentFs = []  # base complex component faces, component face id
        self.componentEs = []  # base complex component edges, component edge id
        self.componentVs = []  # base complex component vertices, component vert id
        self.separated_face_patches = []  # mesh face id
        self.color_num = 8  # for hvs

    def build_face(self):
        print("extracting base complex faces ...")
        self.mesh.clean_visited()
        self.singularity_structure.clean_singular_element_visited()
        for se in self.singularity_structure.singular_Es:
            se.visited = True
            for eid in se.es_link:
                # each face of a singular edge will be assigned to a component
                edge = self.mesh.Edges[eid]
                for fid in edge.neighboring_Fids:
                    if self.mesh.Faces[fid].visited:
                        continue
                    pathFids = self.mesh.trace_patch_faces(fid)
                    self.separated_face_patches.append(pathFids)  # not be used

        for fid, face in self.mesh.Faces.items():
            if face.is_super:
                continue
            if face.visited or face.is_boundary:
                self.Fids.append(fid)
        print("done ...")

    # https://github.com/Cotrik/CotrikMesh/blob/d890c9e596fa3242342b8a1b7fc7b9ce194f03a7/libcotrik/src/BaseComplex.cpp#L1420
    def build_edge(self):
        print("extracting base complex edges ...")
        if not self.Fids:
            print("No Fids in base complex ...")
            return
        self.mesh.clean_visited()
        # self.singularity_structure.cleanSingularElementVisited()
        tmp_eids = set()
        # mark all edges of base complex faces
        # only keep the eids that all neighboring face in Fids
        # Fids means be used by a base complex patch face or is_boundary.
        for fid in self.Fids:
            face = self.mesh.Faces[fid]
            for eid in face.Eids:
                tmp_eids.add(eid)
        remove_eids = set()
        for eid in tmp_eids:
            edge = self.mesh.Edges[eid]
            for nfid in edge.neighboring_Fids:
                if nfid not in self.Fids:
                    remove_eids.add(eid)
        # remove
        for e in remove_eids:
            tmp_eids.remove(e)
        self.Eids = list(tmp_eids)
        print("done ...")

    # https://github.com/Cotrik/CotrikMesh/blob/d890c9e596fa3242342b8a1b7fc7b9ce194f03a7/libcotrik/src/BaseComplex.cpp#L1026
    def build_vertex(self):
        print("extracting base complex vertex ...")
        if not self.Fids or not self.Eids:
            print("No Fids or Eids in base complex ...")
            return
        tmp_vids = set([])
        # find all vertices that used in at least one base complex Fid
        for fid in self.Fids:
            face = self.mesh.Faces[fid]
            for vid in face.Vids:
                tmp_vids.add(vid)
        # keep vid that all neighboring eid are in Eids
        remove_vids = set()
        for vid in tmp_vids:
            vert = self.mesh.Vertices[vid]
            for e in vert.neighboring_Eids:
                if e not in self.Eids:
                    remove_vids.add(vid)
        # remove
        for v in remove_vids:
            tmp_vids.remove(v)
        self.Vids = list(tmp_vids)
        print("done ...")

    def region_grow_edge(self, start_eid):
        eid_stack = [start_eid]
        tmp_eids = set()
        while eid_stack:
            cur_eid = eid_stack.pop()
            cur_edge = self.mesh.Edges[cur_eid]
            if cur_edge.visited:
                continue
            cur_edge.visited = True
            tmp_eids.add(cur_eid)
            for vid in cur_edge.Vids:
                if vid in self.Vids:
                    continue
                # get next edge
                next_eids = self.mesh.get_next_edge(vid, cur_eid)
                eid_stack.extend(next_eids)
        return list(tmp_eids)

    # in this function plz do not reset the visited flag of mesh
    # the visited flag will be used in the parent function
    def trace_vert_for_vid_link_and_eid_link_component_edge(self, start_vid, component_edge):
        # goal is to get a ordered edge id link
        # self.mesh.clean_visited()
        vid_visited = [False] * len(self.mesh.Vertices)
        eid_visited = [False] * len(self.mesh.Edges)
        vids_link = []
        eids_link = []
        vids_stack = [start_vid]
        while vids_stack:
            cur_vid = vids_stack.pop()
            if vid_visited[cur_vid]:
                continue
            vid_visited[cur_vid] = True
            vids_link.append(cur_vid)
            for neid in self.mesh.Vertices[cur_vid].neighboring_Eids:
                if neid not in component_edge.es_link:
                    continue
                cur_edge = self.mesh.Edges[neid]
                if eid_visited[neid]:
                    continue
                eid_visited[neid] = True
                eids_link.append(neid) # add current eid to the link
                # next vid
                for nvid in cur_edge.Vids:
                    if nvid == cur_vid:
                        continue
                    vids_stack.append(nvid)
        return vids_link, eids_link

    # https://github.com/Cotrik/CotrikMesh/blob/d890c9e596fa3242342b8a1b7fc7b9ce194f03a7/libcotrik/src/BaseComplex.cpp#L1559
    def build_component_edge_link(self, component_edge):
        all_evids = []
        for eid in component_edge.es_link:
            eVids = self.mesh.Edges[eid].Vids
            all_evids.extend(eVids)
        unique_evids = set(all_evids)
        start_vid = all_evids[0]
        if len(unique_evids) != len(component_edge.es_link):
            # not circular
            # however, if there contain sub loop
            # the condition is not correct
            # find start vid
            start_flag = [False] * len(self.mesh.Vertices)
            for v in all_evids:
                start_flag[v] = not start_flag[v]
            # first True vid is the start vid
            for v in all_evids:
                if start_flag[v]:
                    start_vid = v
                    break
        component_edge.vs_link = list(unique_evids)
        # start vertex
        start_vert = self.mesh.Vertices[start_vid]
        start_eid = -1
        for neid in start_vert.neighboring_Eids:
            if neid in component_edge.es_link:
                start_eid = neid
                break
        start_edge = self.mesh.Edges[start_eid]
        # travel for component edge's vid link and eid link
        vids_link, eids_link = self.trace_vert_for_vid_link_and_eid_link_component_edge(start_vid, component_edge)
        component_edge.vs_link = vids_link
        component_edge.es_link = eids_link

    def build_component_edges(self):
        print("building base complex component edge ...")
        self.mesh.clean_visited()
        for eid in self.Eids:
            edge = self.mesh.Edges[eid]
            if edge.visited:
                continue
            # region
            # when build base complex,
            # a mesh edge only could be traced once
            component_eid_link = self.region_grow_edge(eid) # already sorted eid list
            ceid = len(self.componentEs)
            for e in component_eid_link:
                self.mesh.Edges[e].component_ids.add(ceid)
            ce = ComponentEdge(ceid, component_eid_link, edge.is_boundary)
            self.build_component_edge_link(ce)
            self.componentEs.append(ce)
        print(f'done ... {len(self.componentEs)} Component Edges')

    #                 ---------
    #                 |        |
    #                 | * < -- | ----  currentFace
    #                 |        |
    # currentEdge --> ---------
    #                 |        |
    #                 | * < -- | ----  nextFace
    #                 |        |
    #                 ---------
    def region_grow_face(self, start_fid):
        fid_stack = [start_fid]
        tmp_fids = set(fid_stack)
        while fid_stack:
            cur_fid = fid_stack.pop()
            cur_face = self.mesh.Faces[cur_fid]
            if cur_face.visited:
                continue
            cur_face.visited = True
            tmp_fids.add(cur_fid)
            for eid in cur_face.Eids:
                if eid in self.Eids:
                    continue
                # edge = self.mesh.Edges[eid]
                oppsite_face = self.mesh.get_opposite_face_on_any_type_eid(cur_fid, eid)
                fid_stack.extend(oppsite_face)
        return list(tmp_fids)

    def build_component_faces(self):
        print("building base complex component face ...")
        self.mesh.clean_visited()
        # self.singularity_structure.cleanSingularElementVisited()
        for fid in self.Fids:
            face = self.mesh.Faces[fid]
            if face.visited:
                continue
            region_fids = self.region_grow_face(fid)
            componentId = len(self.componentFs)
            for rfid in region_fids:
                self.mesh.Faces[rfid].component_ids.add(componentId)
            cFace = ComponentFace(componentId, region_fids, face.is_boundary)
            self.componentFs.append(cFace)
        print(f'done ... {len(self.componentFs)} Component Faces')

    def region_grow_cell(self, start_cid):
        cid_stack = [start_cid]
        regionCids = set()
        while cid_stack:
            cur_cid = cid_stack.pop()
            cur_c = self.mesh.Cells[cur_cid]
            if cur_c.visited:
                continue
            regionCids.add(cur_cid) # add to region cids
            cur_c.visited = True
            fids = cur_c.Fids
            for fid in fids:
                if fid in self.Fids:
                    continue
                ncids = self.mesh.Faces[fid].neighboring_Cids
                # manifold
                # or non manifold ncids > 2
                for ncid in ncids:
                    if ncid == cur_cid:
                        continue
                    cid_stack.append(ncid)
        return list(regionCids)

    def build_component_cells(self):
        print("building base complex component cell ...")
        self.mesh.clean_visited()

        for cid, cell in self.mesh.Cells.items():
            if cell.is_super:
                continue
            if cell.visited:
                continue
            regionCells = self.region_grow_cell(cid)
            componentCid = len(self.componentCs)
            for rcid in regionCells:
                self.mesh.Cells[rcid].component_ids.add(componentCid)
            # component cell
            cc = ComponentCell(componentCid, regionCells)
            self.componentCs.append(cc)
        print(f'done ... {len(self.componentCs)} Component Cells')

    def build_component_vertices(self):
        print("building base complex component vertex ...")
        for i, vid in enumerate(self.Vids):
            self.mesh.Vertices[vid].component_ids.add(i)
            cv = ComponentVert(i, mesh_vid=vid, xyz=self.mesh.Vertices[vid].xyz())
            cv.is_boundary = self.mesh.Vertices[vid].is_boundary
            self.componentVs.append(cv)
        print(f'done ... {len(self.componentVs)} Component Vertices')

    def build_components(self):
        self.build_component_cells()
        # component color
        self.build_component_faces()
        self.build_component_edges()
        self.build_component_vertices()

    # ---------------------------------------------------------------------
    # build component connectivities
    # https://github.dev/Cotrik/CotrikMesh/blob/master/libcotrik/src/BaseComplex.cpp
    # ---------------------------------------------------------------------
    # get the fids that belongs to component cell and component face
    def get_fids_of_component_cell(self, component_cell):
        constituent_component_fids = set()
        boundary_mesh_ids = set()
        for cid in component_cell.cids_patch:
            cfids = self.mesh.Cells[cid].Fids
            for f in cfids:
                f_component_fid = self.mesh.Faces[f].component_ids
                constituent_component_fids.update(f_component_fid)
                if self.mesh.is_face_at_boundary_of_component_cell(f):
                    boundary_mesh_ids.add(f)
        # find out the real constituent fids
        constituent_component_fids = self.find_constituent_fids(constituent_component_fids, boundary_mesh_ids)
        # varify correctness
        all_candidate_mesh_fids = set()
        for fid in constituent_component_fids:
            face = self.componentFs[fid]
            all_candidate_mesh_fids.update(face.fids_patch)
        if all_candidate_mesh_fids != boundary_mesh_ids:
            print(f'component cell {component_cell.element_id} surface is not closed.')
        if len(constituent_component_fids) == 0:
            print(f'component cell {component_cell.element_id} does not have boundary faces.')
        return list(constituent_component_fids)

    def find_constituent_fids(self, candidate_component_fids, bnd_mesh_fids):
        candidate_fids = candidate_component_fids
        # fid is component face id
        for fid in candidate_component_fids:
            face = self.componentFs[fid]
            out_scope_fids = set(face.fids_patch) - bnd_mesh_fids
            if out_scope_fids:
                continue
            new_candidate_fids = []
            for cf in candidate_fids:
                bigger_fid = self.get_bigger_intersecting_component_face(fid, cf)
                if bigger_fid:
                    new_candidate_fids.append(bigger_fid)
                else:
                    new_candidate_fids.append(cf)
            candidate_fids = new_candidate_fids
        return candidate_fids

    # compair size of two overlapping component faces
    # if two component faces does not have any common mesh face,
    #   return None
    # else
    #   return larger component face id
    def get_bigger_intersecting_component_face(self, fid_0, fid_1):
        if fid_0 == fid_1:
            return None
        face_patch_0 = set(self.componentFs[fid_0].fids_patch)
        face_patch_1 = set(self.componentFs[fid_1].fids_patch)
        if face_patch_0.intersection(face_patch_1):
            if __debug__ and (len(face_patch_0) == len(face_patch_1)):
                print(f'component face {fid_0} contain same mesh face patch with component face {fid_1}.')
            return fid_0 if len(face_patch_0) >= len(face_patch_1) else fid_1
        return None

    # return component vids
    def get_vids_of_component_face(self, component_face):
        faceVids = set()
        for ce in component_face.Eids:
            ids = self.componentEs[ce].Vids
            if len(ids) != 0:
                faceVids.update(ids)
        if __debug__ and (len(faceVids) != 4):
            print("component face {} contain {} vertices.".format(component_face.element_id, len(faceVids)))
        if len(faceVids) != 4:
            return faceVids
        # if the face is quad, continue
        # build face vids order
        eid = component_face.Eids[0]
        v1 = self.componentEs[eid].Vids[0]
        v2 = self.componentEs[eid].Vids[1]
        v3 = -1
        # v3
        for ce in component_face.Eids[1:]:
            ids = self.componentEs[ce].Vids
            if ids[0] == v2:
                v3 = ids[1]
            if ids[1] == v2:
                v3 = ids[0]
        v4 = -1
        for ce in component_face.Eids[1:]:
            ids = self.componentEs[ce].Vids
            if ids[0] == v3 and ids[1] != v2:
                v4 = ids[1]
            if ids[1] == v3 and ids[0] != v2:
                v4 = ids[0]
        return [v1, v2, v3, v4]

    def get_vids_of_hexa_component_cell(self, component_cell):
        top_face_vids = self.componentFs[component_cell.Fids[0]].Vids
        bottom_face_vids = []
        final_vids = set(top_face_vids)
        # hexa
        for fid in component_cell.Fids:
            vids = self.componentFs[fid].Vids
            top_bottom_face_vids = set(top_face_vids)
            top_bottom_face_vids.update(vids)
            if len(top_bottom_face_vids) == 8:
                bottom_face_vids = vids
                final_vids.update(vids)
                break
        if not bottom_face_vids:  # could not find bottom face vids
            print("circle cell ... not vid is returned")
            return []
        return final_vids

    def get_vids_of_non_hexa_component_cell(self, component_cell):
        # each non hexa cell are a cell component
        final_vids = set()
        for fid in component_cell.Fids:
            vids = self.componentFs[fid].Vids
            final_vids.update(vids)
        return final_vids

    # return component vids
    def get_vids_of_component_cell(self, component_cell):
        hexa = len(component_cell.Fids) == 6
        if hexa:
            for fid in component_cell.Fids:
                if len(self.componentFs[fid].Vids) != 4:
                    hexa = False
                    break
        if hexa:
            final_vids = self.get_vids_of_hexa_component_cell(component_cell)
        else:
            final_vids = self.get_vids_of_non_hexa_component_cell(component_cell)
        # non hexa
        # in vtk file, face need to be order by a special order
        # in blender, a cell is rendered face by face, so, vid order does not matter.
        return final_vids

    def get_neighbor_cids_for_component_cell(self, component_cell):
        neighbor_cids = set()
        for fid in component_cell.Fids:
            cf_ncids = self.componentFs[fid].neighboring_Cids
            for ncid in cf_ncids:
                if ncid == component_cell.element_id:
                    continue
                neighbor_cids.add(ncid)
        return list(neighbor_cids)

    def assign_color_to_component_cell(self, component_cell):
        neighbor_colors = set()
        for ncid in component_cell.neighboring_Cids:
            cc = self.componentCs[ncid]
            neighbor_colors.add(cc.color)
        if len(neighbor_colors) >= self.color_num - 1:
            self.color_num = len(neighbor_colors) + 2 # extra one for the input cell
        candidate_colors = set()
        for i in range(self.color_num):
            if i == 0:
                continue
            if i not in neighbor_colors:
                candidate_colors.add(i)
        picked_color = random.choice(list(candidate_colors))
        # pickColor = list(candidate_colors)[int(len(candidate_colors) / 2)]
        component_cell.color = picked_color

    def get_eids_of_component_face(self, component_face):
        eids = set()
        all_bnd_mesh_eids = set()
        for fid in component_face.fids_patch:
            (component_eids,
             bnd_mesh_eids) = self.mesh.find_all_bnd_component_eids_of_component_face(fid,
                                                                                      component_face.element_id)
            eids.update(component_eids)
            all_bnd_mesh_eids.update(bnd_mesh_eids)
        # find real bnd component eids
        eids = self.get_constituent_eids(eids, all_bnd_mesh_eids)
        return list(eids)

    # real bnd eids must covert all mesh edges
    def get_constituent_eids(self, candidate_eids, bnd_mesh_eids):
        eids = candidate_eids
        for eid in candidate_eids:
            edge = self.componentEs[eid]
            out_scope_eids = set(edge.es_link) - bnd_mesh_eids
            if out_scope_eids:
                continue
            new_eids = set()
            for ce in eids:
                bigger_eid = self.get_bigger_intersecting_component_eid(eid, ce)
                if bigger_eid:
                    new_eids.add(bigger_eid)
                else:
                    new_eids.add(ce)
            eids = new_eids
            # need to verify 21/aug/2023
        return eids

    def get_bigger_intersecting_component_eid(self, eid_0, eid_1):
        if eid_0 == eid_1:
            return None
        mesh_eid_0 = set(self.componentEs[eid_0].es_link)
        mesh_eid_1 = set(self.componentEs[eid_1].es_link)
        if mesh_eid_0.intersection(mesh_eid_1):
            if len(mesh_eid_0) == len(mesh_eid_1):
                print(f'component edge {eid_0} contain same mesh edge link with component edge {eid_1}.')
            return eid_0 if len(mesh_eid_0) >= len(mesh_eid_1) else eid_1
        return None

    def get_eids_of_component_cell(self, component_cell):
        component_eids = set()
        for cfid in component_cell.Fids:
            cf = self.componentFs[cfid]
            ceids = self.get_eids_of_component_face(cf)
            component_eids.update(ceids)
        return list(component_eids)

    def get_neighbor_fids_of_face_component(self, component_face):
        component_eids = component_face.Eids
        neighbor_component_fids = set()
        for eid in component_eids:
            nfids = self.componentEs[eid].neighboring_Fids
            neighbor_component_fids.update(nfids)
        neighbor_component_fids.remove(component_face.element_id)
        return list(neighbor_component_fids)

    def get_neighbor_cids_of_face_component(self, component_face):
        mesh_fid = component_face.fids_patch[0]
        ncids = self.mesh.Faces[mesh_fid].neighboring_Cids
        # assume mesh is manifold
        neighboring_component_cids = set()
        for ncid in ncids:
            neighboring_component_cids.update(self.mesh.Cells[ncid].component_ids)
        return list(neighboring_component_cids)

    def get_neighbor_cids_of_edge_component(self, component_edge):
        ncids = self.mesh.Edges[component_edge.es_link[0]].neighboring_Cids
        ccids = set()
        for ncid in ncids:
            ccid = self.mesh.Cells[ncid].component_ids
            ccids.update(ccid)
        return list(ccids)

    def get_neighbor_fids_of_edge_component(self, component_edge):
        nfids = self.mesh.Edges[component_edge.es_link[0]].neighboring_Fids
        cfids = set()
        for f in nfids:
            neighbor_component_fids = self.mesh.Faces[f].component_ids
            cfids.update(neighbor_component_fids)
        return list(cfids)

    # vertex should only have one component id
    def get_vids_of_component_edge(self, ce):
        if ce.vs_link[0] == ce.vs_link[-1]:
            return []
        else:
            v0 = self.mesh.Vertices[ce.vs_link[0]].component_ids
            v1 = self.mesh.Vertices[ce.vs_link[-1]].component_ids
            if len(v0) > 1 or len(v1) > 1:
                print("duplication vertex component id ...")
                print(f'v0 : {v0}')
                print(f'v1 : {v1}')
            # https://stackoverflow.com/questions/59825/how-to-retrieve-an-element-from-a-set-without-removing-it
            for e1 in v0:
                break
            v0 = e1
            for e2 in v1:
                break
            v1 = e2
            return [v0, v1]

    def get_component_eids_of_singular_edge(self, se):
        if not se.componentEids_link:
            return se.componentEids_link
        ce_link = set() # component eids link
        for eid in se.es_link:
            ceid = self.mesh.Edges[eid].component_ids
            ce_link.update(ceid)
        se.componentEids_link = ce_link
        return ce_link

    def get_neighbor_fids_of_singular_edge(self, se):
        ce_link = self.get_component_eids_of_singular_edge(se)
        nfids = set() # component face ids
        for ceid in ce_link:
            ce = self.componentEs[ceid]
            tmp_nfids = ce.neighboring_Fids
            nfids.update(tmp_nfids)
        return list(nfids)

    def get_neighbor_cids_of_singular_edge(self, se):
        ce_link = self.get_component_eids_of_singular_edge(se)
        ncids = set() # component cell ids
        for ceid in ce_link:
            ce = self.componentEs[ceid]
            tmp_ncids = ce.neighboring_Cids
            ncids.update(tmp_ncids)
        return list(ncids)

    def build_component_constituent(self):
        print("building component cell [face, edge] constituent ...")
        for cc in self.componentCs:
            # component face id
            cc.Fids = self.get_fids_of_component_cell(cc)
            # component edge id
            cc.Eids = self.get_eids_of_component_cell(cc)
        print("done ...")
        print("building component face [edge] constituent ...")
        for cf in self.componentFs:
            # cf.Eids = self.get_eids_of_component_face(cf)
            cf.update_eids(self.get_eids_of_component_face(cf))
        print("done ...")
        print("building component edge [vertex] constituent ...")
        for ce in self.componentEs:
            ce.Vids = self.get_vids_of_component_edge(ce)
        print("done ...")
        print("building component face [vertex] constituent ...")
        for cf in self.componentFs:
            cf.Vids = self.get_vids_of_component_face(cf)
        print("done ...")
        print("building component Cell [vertex] constituent ...")
        for cc in self.componentCs:
            cc.Vids = self.get_vids_of_component_cell(cc)
        print("done ...")

    def build_component_connectivities(self):
        print("building component cell connectivites ...")
        for cc in self.componentCs:
            cc.Fids = self.get_fids_of_component_cell(cc) # component face id
            cc.Eids = self.get_eids_of_component_cell(cc) # component edge id
        print("done ...")
        print("building component face connectivites ...")
        for cf in self.componentFs:
            cf.Eids = self.get_eids_of_component_face(cf)
            cf.neighboring_Cids = self.get_neighbor_cids_of_face_component(cf)
        print("done ...")
        print("building component edge connectivites ...")
        for ce in self.componentEs:
            ce.Vids = self.get_vids_of_component_edge(ce)
            ce.neighboring_Cids = self.get_neighbor_cids_of_edge_component(ce)
            ce.neighboring_Fids = self.get_neighbor_fids_of_edge_component(ce)
        print("done ...")
        print("building component Face component vids and neighbor component fids ...")
        for cf in self.componentFs:
            cf.Vids = self.get_vids_of_component_face(cf)
            cf.neighboring_Fids = self.get_neighbor_fids_of_face_component(cf)
        print("done ...")
        print("building component Cell component vids, neighbor component cids, and assign color ...")
        for cc in self.componentCs:
            cc.Vids = self.get_vids_of_component_cell(cc)
            cc.neighboring_Cids = self.get_neighbor_cids_for_component_cell(cc)
            self.assign_color_to_component_cell(cc)
        print("done ...")
        print("building component vertex neighbor elements ...")
        for ce in self.componentEs:
            for cvid in ce.Vids:
                vert = self.componentVs[cvid]
                vert.add_neighboring_vids([ce.Vids])
                vert.add_neighboring_eids([ce.element_id])
                vert.add_neighboring_fids(ce.neighboring_Fids)
                vert.add_neighboring_cids(ce.neighboring_Cids)
        print("done ...")
        print("building component edge neighbor edges ...")
        for ce in self.componentEs:
            for cvid in ce.Vids:
                vert = self.componentVs[cvid]
                ce.add_neighboring_eids(vert.neighboring_Eids)
        print("done ...")
        print("building singular edge neighbor elements ...")
        for se in self.singularity_structure.singular_Es:
            se.neighboring_Cids = self.get_neighbor_cids_of_singular_edge(se)
            se.neighboring_Fids = self.get_neighbor_fids_of_singular_edge(se)
            # not yet implement
            # GetNeighborComponentFidsGroups
            # GetNeighborComponentCidsGroups
        print("done ...")
    # ---------------------------------------------------------------------

    # ---------------------------------------------------------------------
    # build singular edge separated patches
    # https://github.dev/Cotrik/CotrikMesh/blob/master/libcotrik/src/BaseComplex.cpp
    # the patches that directly releated to singular edge
    # ---------------------------------------------------------------------
    def build_singular_edge_separated_face_patches(self):
        for se in self.singularity_structure.singular_Es:
            edge = self.mesh.Edges[se.es_link[0]]
            for nfid in edge.neighboring_Fids:
                for pid, patch in enumerate(self.separated_face_patches):
                    if nfid in patch:
                        se.separatedFacePatchIds.add(pid)

    def build_singular_edge_separated_component_face_pathces(self):
        pass

    def build_singular_edge_separated_pathces(self):
        pass

    # ---------------------------------------------------------------------
    def build(self):
        self.build_face() # build first
        self.build_edge() # edge is build after face
        self.build_vertex() # using Fids and Eids
        self.build_components()
        self.build_component_constituent()
        # connectivity is build in mesh structure
        # self.build_component_connectivities()

    # ----------------------------------------------------------------------
    # support function
    # ----------------------------------------------------------------------
    def detect_same_vid_faces(self):
        for fid, face in self.mesh.Faces.items():
            for fid_1, face_1 in self.mesh.Faces.items():
                if fid == fid_1:
                    continue
                if set(face.Vids) == set(face_1.Vids):
                    print("fid {} is same with face {}".format(fid, fid_1))

    # ----------------------------------------------------------------------
    # -- detect same component elements
    # ----------------------------------------------------------------------
    def is_same_component_vert(self, component_vert_0, component_vert_1):
        component_vid_0 = component_vert_0.element_id
        component_vid_1 = component_vert_1.element_id
        if component_vid_0 == component_vid_1:
            print("two input vids are same ...")
            return True, True
        vert_0 = component_vert_0
        vert_1 = component_vert_1
        # same xyz
        same_xyz = False
        if vert_0.xyz == vert_1.xyz:
            same_xyz = True
        # same mesh vid
        same_m_vid = False
        if vert_0.mesh_vid == vert_1.mesh_vid:
            same_m_vid = True
        if same_xyz or same_m_vid:
            print("vert {} has same xyz with vert {} : {}".format(component_vid_0, component_vid_1, same_xyz))
            print("vert {} has same mesh id with vert {} : {}".format(component_vid_0, component_vid_1, same_m_vid))
        return same_xyz, same_m_vid

    # id, vids, es_link, vs_link
    def is_same_component_edge(self, edge_0, edge_1):
        component_eid_0 = edge_0.element_id
        component_eid_1 = edge_1.element_id
        if component_eid_0 == component_eid_1:
            print("two input eids are same ...")
            return True, True, True

        # same vids
        same_vids = False
        if set(edge_0.Vids) == set(edge_1.Vids):
            same_vids = True
        same_es_link = False
        if set(edge_0.es_link) == set(edge_1.es_link):
            same_es_link = True
        same_vs_link = False
        if set(edge_0.vs_link) == set(edge_1.vs_link):
            same_vs_link = True
        if same_vids or same_es_link or same_vs_link:
            print("edge {} and edge {}".format(component_eid_0, component_eid_1))
            print("same vids : {}".format(same_vids))
            print("same es link : {}".format(same_es_link))
            print("same vs link : {}".format(same_vs_link))
        return same_vids, same_es_link, same_vs_link

    # id, vids, eids, fid_patch
    def is_same_component_face(self, face_0, face_1):
        component_fid_0 = face_0.element_id
        component_fid_1 = face_1.element_id
        if component_fid_0 == component_fid_1:
            print("two input fids are same ...")
            return True, True, True
        same_vids = False
        if set(face_0.Vids) == set(face_1.Vids):
            same_vids = True
        same_eids = False
        if set(face_0.Eids) == set(face_1.Eids):
            same_eids = True
        same_fid_patch = False
        if set(face_0.fids_patch) == set(face_1.fids_patch):
            same_fid_patch = True
        if same_vids or same_eids or same_fid_patch:
            print("face {} and face {}".format(component_fid_0, component_fid_1))
            print("same vids : {}".format(same_vids))
            print("same eids : {}".format(same_eids))
            print("same fid patch : {}".format(same_fid_patch))
        return same_vids, same_eids, same_fid_patch

    # id, vid, eid, fid, cid_patch
    def is_same_component_cell(self, cell_0, cell_1):
        cid_0 = cell_0.element_id
        cid_1 = cell_1.element_id
        if cid_0 == cid_1:
            print("two input cell are same ...")
            return True, True, True, True
        same_vids = False
        if set(cell_0.Vids) == set(cell_1.Vids):
            same_vids = True
        same_eids = False
        if set(cell_0.Eids) == set(cell_1.Eids):
            same_eids = True
        same_fids = False
        if set(cell_0.Fids) == set(cell_1.Fids):
            same_fids = True
        same_cid_patch = False
        if set(cell_0.cids_patch) == set(cell_1.cids_patch):
            same_cid_patch = True
        if same_vids or same_eids or same_fids or same_cid_patch:
            print("cell {} and cell {}".format(cid_0, cid_1))
            print("same vids : {}".format(same_vids))
            print("same eids : {}".format(same_eids))
            print("same fids : {}".format(same_fids))
            print("same cid patch : {}".format(same_cid_patch))
        return same_vids, same_eids, same_fids, same_cid_patch

    def detect_duplicate_component_vert(self):
        visited = []
        print("detecting duplicate component vert -----")
        for vert in self.componentVs:
            visited.append(vert.element_id)
            for vert_1 in self.componentVs:
                if vert_1.element_id in visited:
                    continue
                if vert.element_id == vert_1.element_id:
                    continue
                self.is_same_component_vert(vert, vert_1)
        print("detecting duplicate component vert ----- done")

    def detect_duplicate_component_edge(self):
        visited = []
        same_eids = []
        print("detecting duplicate component edge -----")
        for edge_0 in self.componentEs:
            visited.append(edge_0.element_id)
            if edge_0.is_super:
                continue
            for edge_1 in self.componentEs:
                if edge_1.element_id in visited:
                    continue
                if edge_1.is_super:
                    continue
                if edge_0.element_id == edge_1.element_id:
                    continue
                v, es, vs = self.is_same_component_edge(edge_0, edge_1)
                if v or es or vs:
                    tmp = [edge_0.element_id, edge_1.element_id]
                    tmp.sort()
                    if tmp not in same_eids:
                        same_eids.append(tmp)
        print("detecting duplicate component edge ----- done")
        return same_eids

    def detect_duplicate_component_face(self):
        visited = []
        same_fids = []
        print("detecting duplicate component face -----")
        for face_0 in self.componentFs:
            visited.append(face_0.element_id)
            if face_0.is_super:
                continue
            for face_1 in self.componentFs:
                if face_1.element_id in visited:
                    continue
                if face_1.is_super:
                    continue
                if face_0.element_id == face_1.element_id:
                    continue
                v, e, f = self.is_same_component_face(face_0, face_1)
                if v or e or f:
                    tmp = [face_0.element_id, face_1.element_id]
                    tmp.sort()
                    if tmp not in same_fids:
                        same_fids.append(tmp)
        print("detecting duplicate component face ----- done")
        return same_fids

    def detect_duplicate_component_cell(self):
        visited = []
        print("detecting duplicate component cell -----")
        for cell_0 in self.componentCs:
            visited.append(cell_0.element_id)
            if cell_0.is_super:
                continue
            for cell_1 in self.componentCs:
                if cell_1.element_id in visited:
                    continue
                if cell_1.is_super:
                    continue
                if cell_0.element_id == cell_1.element_id:
                    continue
                self.is_same_component_cell(cell_0, cell_1)
        print("detecting duplicate component cell ----- done")

    def get_vertex_list(self):
        vertex_list = []
        for vid in self.componentVs:
            coord = vid.xyz()
            vertex_list.append(coord)
        return vertex_list
