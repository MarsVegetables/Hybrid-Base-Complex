from blender_viz.blender_viz_basic_op import *
from mesh_structure.mesh_common import check_common_elements


class BlenderMeshVisualizer:
    def __init__(self, mesh):
        self.mesh = mesh
        # self.mesh = SuperHybridMesh(None, None)
        self.vertices = self.mesh.get_vertex_list()
        self.cleanEdges = self.mesh.get_edge_list(skipRemoved=True)
        self.cleanFaces = self.mesh.get_face_list(skipRemoved=True)
        self.default_r = self.mesh.get_average_edge_length() / 15
        self.edge_template = None
        self.sphere_tmp = None

    def show_transparent_boundary(self,
                                  rgba=(0.9, 0.9, 0.9, 0.15),
                                  collection_name="transparent boundary",
                                  cleaned_bnd=True):
        if cleaned_bnd:
            boundary = self.mesh.get_manifold_boundary_fids()
        else:
            boundary = self.mesh.get_boundary_fids()
        self.show_face(fids=boundary,
                       collection_name=collection_name,
                       rgba=rgba)

    def show_wireframe(self, collection_name="mesh wireframe"):
        create_wireframe_geo_node(curve_radius=self.default_r * 0.2)
        obj_name = make_hybrid_mesh_wireframe(self.vertices, self.cleanEdges,
                                              self.cleanFaces, collection_name=collection_name)
        return obj_name

    def is_multi_level_structure(self):
        return self.mesh.structure_type == "multi level hybrid"

    # what is the edge format
    def show_edges(self, eids, edge_names=None,
                   radii=0, collection_name='edge',
                   rgba=(1, 0, 0, 1), show_ground_mesh_edge=True,
                   colored_vids=None, adjust_alpha=False):
        edge_radii = self.default_r * 0.5
        if radii > 0:
            edge_radii = radii
        if not self.edge_template:
            self.edge_template = create_cylinder_template()
        if not edge_names:
            use_edge_names = [f'edge {i}' for i in eids]
        else:
            use_edge_names = edge_names
        if self.is_multi_level_structure():
            (current_mesh,
             edge_vids,
             current_edge_names) = self.mesh.covert_to_background_edge_vids(eids, use_edge_names,
                                                                            show_ground_mesh_edge)
        else:
            current_mesh = self.mesh
            edge_vids = self.mesh.eids_to_edge_vids(eids)
            current_edge_names = use_edge_names
        # highlight the edge connected to the colored vertices
        alpha_01_edge_vids = []
        alpha_01_edge_names = []
        alpha_1_edge_vids = []
        alpha_1_edge_names = []
        if colored_vids:
            _, mesh_vids = self.mesh.covert_to_background_vids(colored_vids)
            for i, ev in enumerate(edge_vids):
                if check_common_elements(ev, mesh_vids):
                    alpha_1_edge_vids.append(ev)
                    alpha_1_edge_names.append(current_edge_names[i])
                    continue
                # other edge are alpha is 0.1
                alpha_01_edge_names.append(current_edge_names[i])
                alpha_01_edge_vids.append(ev)
        elif adjust_alpha:
            alpha_01_edge_vids = edge_vids
            alpha_01_edge_names = current_edge_names
        else:
            alpha_1_edge_vids = edge_vids
            alpha_1_edge_names = current_edge_names

        if alpha_01_edge_vids:
            tmp_rgba = [x for x in rgba]
            tmp_rgba[-1] = 0.2 # in paper
            # tmp_rgba[-1] = 0.1  # in auto render
            # else, all edge alpha are 1
            make_edges(current_mesh.get_vertex_list(), alpha_01_edge_vids, r=edge_radii,
                       name=collection_name, rgba=tmp_rgba, edge_names=alpha_01_edge_names,
                       edge_template=self.edge_template)
        if alpha_1_edge_vids:
            # else, all edge alpha are 1
            make_edges(current_mesh.get_vertex_list(), alpha_1_edge_vids, r=edge_radii,
                       name=collection_name, rgba=rgba, edge_names=alpha_1_edge_names,
                       edge_template=self.edge_template)

    # def show_edges_poly_line(self, eids, edge_names=None,
    #                          radii=0, collection_name='edge',
    #                          rgba=(1, 0, 0, 1), show_ground_mesh_edge=True):
    #     edge_radii = self.default_r * 0.5
    #     if radii > 0:
    #         edge_radii = radii
    #     if not self.edge_template:
    #         self.edge_template = create_cylinder_template()
    #     if not edge_names:
    #         use_edge_names = [f'edge {i}' for i in eids]
    #     else:
    #         use_edge_names = edge_names
    #     if self.is_multi_level_structure():
    #         (current_mesh,
    #          edge_vids,
    #          current_edge_names) = self.mesh.covert_to_background_edge_vids(eids, use_edge_names,
    #                                                                         show_ground_mesh_edge)
    #     else:
    #         current_mesh = self.mesh
    #         edge_vids = self.mesh.eids_to_edge_vids(eids)
    #         current_edge_names = use_edge_names
    #     make_edges(current_mesh.get_vertex_list(), edge_vids, r=edge_radii,
    #                name=collection_name, rgba=rgba, edge_names=current_edge_names,
    #                edge_template=self.edge_template)

    # all input faces will be shown in same mesh
    # so only have one face name of the input faces
    def show_face(self, fids, face_name='face',
                  collection_name='faces', rgba=(1, 1, 1, 1),
                  show_ground_mesh_face=True):
        if self.is_multi_level_structure():
            current_mesh, face_vids = (
                self.mesh.covert_fids_to_background_face_vids(fids=fids,
                                                              show_ground_mesh_face=show_ground_mesh_face)
            )
            current_face_name = face_name
        else:
            current_mesh = self.mesh
            face_vids = self.mesh.fids_to_face_vids(fids)
            current_face_name = face_name
        visualize_faces(vertices=current_mesh.get_vertex_list(),
                        face_vid_groups=face_vids,
                        face_name=current_face_name,
                        rgba=rgba,
                        collection_name=collection_name
                        )

    def show_cells(self, cids, face_name='cells',
                   collection_name="cells", show_ground_mesh_cell=True,
                   rgba=None, alpha=1):
        if self.is_multi_level_structure():
            if rgba and rgba[-1] != 1:
                cell_fids = self.mesh.cids_to_fids(cids, only_boundary=True)
                self.show_face(fids=cell_fids, face_name=f'{face_name}_transparent',
                               collection_name=collection_name + "_transparent", rgba=rgba,
                               show_ground_mesh_face=show_ground_mesh_cell)
                return
            for cid in cids:
                cell_fids = self.mesh.Cells[cid].Fids
                if rgba:
                    cell_rgba = rgba
                else:
                    cell_rgba = color_id_to_rgba(color_id=self.mesh.Cells[cid].color,
                                                 total_color=self.mesh.color_num, alpha=alpha)
                self.show_face(fids=cell_fids, face_name=f'{face_name}-cell {cid}'
                               , collection_name=collection_name, rgba=cell_rgba,
                               show_ground_mesh_face=show_ground_mesh_cell)
        else:
            current_mesh = self.mesh
            cell_face_vids = self.mesh.cids_to_face_vids(cids)
            visualize_faces(vertices=current_mesh.get_vertex_list(),
                            face_vid_groups=cell_face_vids,
                            face_name=face_name,
                            rgba=rgba,
                            collection_name=collection_name
                            )

    # vert name = [v1, v2, v3, ...]
    def show_vertices(self, vids, vert_name=None, collection_name="vertices",
                      rgba=(1, 0, 0, 1), radii=0, radii_weight=None):
        using_vert_name = vert_name
        if vert_name is None:
            using_vert_name = [f'vert {i}' for i in vids]
        if not self.sphere_tmp:
            if radii == 0:
                self.sphere_tmp = create_sphere_tmp(self.default_r * 1.5)
            else:
                self.sphere_tmp = create_sphere_tmp(radii)
        using_r = radii
        if radii == 0:
            using_r = self.default_r
        for i, vid in enumerate(vids):
            xyz_co = self.mesh.Vertices[vid].xyz()
            # print(xyz_co)
            x = 1
            if radii_weight:
                x = radii_weight[vid]
            make_sphere(xyz_co, radius=x * using_r, sphere_name=using_vert_name[i],
                        collection_name=collection_name, rgba=rgba, size_weight=x,
                        sphere_template=self.sphere_tmp)




