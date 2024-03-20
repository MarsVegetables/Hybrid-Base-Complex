from mesh_structure.basic_mesh import BasicMesh


# no super element in this structure
class MultiLevelHybridMesh(BasicMesh):
    def __init__(self, file_path, mesh_type, mesh_level=0):
        BasicMesh.__init__(self, file_path, mesh_type)
        self.structure_type = "multi level hybrid"
        # mesh_level is 0 means background mesh
        self.mesh_level = mesh_level
        # root mesh does not have an input mesh
        self.input_mesh = None

    def build(self):
        # possible situation
        # two end is same but the link path is not same
        self.build_vertex_connectivity()
        self.build_element_connectivity()
        # find out all neighbor relationship

        # build partial neighbor relationship

        # based on face connectivity to mark is_boundary elements
        self.mark_boundary_elements()
        #
        self.mark_all_cell_type()
        #
        self.build_opposite_face_relationship()
        self.build_parallel_face_relationship()
        self.build_parallel_edge_relationship()
        #
        self.assign_color_to_cells()

    def covert_to_background_vids(self, vids, show_ground_mesh_edge=True):
        current_vids = vids
        current_mesh = self
        if show_ground_mesh_edge:
            # map to ground mesh, the iteration should be 0
            while current_mesh.mesh_level > 0:
                new_vids = set()
                for i, v in enumerate(current_vids):
                    new_vids.add(current_mesh.Vertices[v].mesh_vid)

                # assign lower level vids to current vids
                current_vids = new_vids
                # move to lower level mesh
                current_mesh = current_mesh.input_mesh
        return current_mesh, list(current_vids)

    def covert_to_background_edge_vids(self, eids, edge_names=['edge'], show_ground_mesh_edge=True):
        assert len(eids) == len(edge_names)
        current_eids = eids
        current_edge_names = edge_names
        current_mesh = self
        if show_ground_mesh_edge:
            # map to ground mesh, the iteration should be 0
            while current_mesh.mesh_level > 0:
                new_edge_names = []
                new_eids = []
                for i, e in enumerate(current_eids):
                    new_eids.extend(current_mesh.Edges[e].es_link)
                    for j in current_mesh.Edges[e].es_link:
                        new_edge_names.append(current_edge_names[i])
                # assign lower level eids to current eids
                # note : [[eid 1, ...], [eid 2, ...], ...]
                current_eids = new_eids
                current_edge_names = new_edge_names
                # move to lower level mesh
                current_mesh = current_mesh.input_mesh
        # covert eid to edge vids [v1, v2]
        edge_vids = current_mesh.eids_to_edge_vids(current_eids)
        return current_mesh, edge_vids, current_edge_names

    def covert_to_background_eids(self, eids, edge_names=['edge'], show_ground_mesh_edge=True):
        assert len(eids) == len(edge_names)
        current_eids = eids
        current_edge_names = edge_names
        current_mesh = self
        if show_ground_mesh_edge:
            # map to ground mesh, the iteration should be 0
            while current_mesh.mesh_level > 0:
                new_edge_names = []
                new_eids = []
                for i, e in enumerate(current_eids):
                    new_eids.extend(current_mesh.Edges[e].es_link)
                    for j in current_mesh.Edges[e].es_link:
                        new_edge_names.append(current_edge_names[i])
                # assign lower level eids to current eids
                # note : [[eid 1, ...], [eid 2, ...], ...]
                current_eids = new_eids
                current_edge_names = new_edge_names
                # move to lower level mesh
                current_mesh = current_mesh.input_mesh
        return current_mesh, current_eids, current_edge_names

    def covert_fids_to_background_face_vids(self, fids, face_names=['face'],
                                            show_ground_mesh_face=True):
        # assert len(fids) == len(face_names)
        current_ids = fids
        # current_names = face_names
        current_mesh = self
        if show_ground_mesh_face:
            # map to ground mesh, the iteration should be 0
            while current_mesh.mesh_level > 0:
                # new_names = []
                new_ids = []
                for i, e in enumerate(current_ids):
                    new_ids.extend(current_mesh.Faces[e].fids_patch)
                    # new_names.append(current_names[i])
                # assign lower level ids to current ids
                # note : [[element 1 - lower id 1, ...],
                #         [element 2 - lower id 1, ...],
                #         ...]
                current_ids = new_ids
                # current_names = new_names
                # move to lower level mesh
                current_mesh = current_mesh.input_mesh
        # covert element id to element vids [v1, v2, ...]
        element_vids = current_mesh.fids_to_face_vids(current_ids)
        # return current_mesh, element_vids, current_names
        return current_mesh, element_vids

    # -------------------------------------------------------------------
    # build all neighbor cell connectivity for face
    # must be component elements
    # -------------------------------------------------------------------
    # def find_all_neighbor_cells_of_component_face(self):
    #     if self.mesh_level == 0:
    #         return
    #     for fid, component_face in self.Faces.items():
    #         if not component_face.is_component():
    #             continue
    #         # for mesh_idxxx
    #         mesh_fid = component_face.fids_patch[0]
    #         ncids = self.input_mesh.Faces[mesh_fid].neighboring_Cids
    #         # assume mesh is manifold
    #         neighboring_component_cids = set()
    #         for ncid in ncids:
    #             neighboring_component_cids.add(self.input_mesh.Cells[ncid].component_id)
    #         return list(neighboring_component_cids)

    # -------------------------------------------------------------------
    # find component boundary element
    # must be component elements
    # -------------------------------------------------------------------
    def is_face_at_boundary_of_component_cell(self, fid):
        face = self.Faces[fid]
        if face.is_boundary:
            return True
        component_cids = set()
        neighboring_cids = face.neighboring_Cids_a if face.neighboring_Cids_a else face.neighboring_Cids
        for cid in neighboring_cids:
            cell = self.Cells[cid]
            component_cids.update(cell.component_ids)
        if len(component_cids) > 1:
            return True
        return False

    # for get_eids_of_component_face
    def find_all_bnd_component_eids_of_component_face(self, fid, component_fid):
        # may include T edge
        bnd_component_eids = set()
        bnd_mesh_eids = set()
        face_eids = self.Faces[fid].Eids
        for eid in face_eids:
            # count :
            # counts the number of neighboring faces that belongs to input component fid
            # count > 1 must be an inner edge
            count = 0
            neighbor_fids = self.Edges[eid].neighboring_Fids
            for nfid in neighbor_fids:
                if component_fid in self.Faces[nfid].component_ids:
                    count += 1
            if count == 1:
                bnd_mesh_eids.add(eid)
                bnd_component_eids.update(self.Edges[eid].component_ids)
        return bnd_component_eids, bnd_mesh_eids

    # -------------------------------------------------------------------
    # construct all relationship
    # must be component elements
    # -------------------------------------------------------------------
    def build_all_connectivity(self):
        if self.mesh_level == 0:
            return
        pass

    # -------------------------------------------------------------------
    # build partial relationship
    # must be component elements
    # -------------------------------------------------------------------
    def build_partial_connectivity(self):
        if self.mesh_level == 0:
            return
        pass

    # vertex does not have any partial neighboring vertex
    # but could have partial neighboring edge, face, and cell
    def build_vertex_partial_connectivity(self):
        pass

    def build_edge_partial_connectivity(self):
        pass

    def build_face_partial_connectivity(self):
        pass

    # if two cells share any cell
    # if two cells share partial faces
    # if two cell share partial edges
    # if two cell share partial vertex
    def build_cell_partial_connectivity(self):
        for fid, face in self.Faces.items():
            pass

    def set_edges_resolution(self):
        if self.mesh_level == 0:
            return
        for eid, edge in self.Edges:
            edge.resolution = len(edge.es_link)
