from blender_viz.mesh_visualizer import BlenderMeshVisualizer
from blender_viz.blender_viz_basic_op import make_edges
# from hybrid_base_complex.hybrid_singularity_structure import HybridSingularityGraph
# from blender_viz.blender_viz_basic_op import make_poly_line


class SingularityGraphVisualizer(BlenderMeshVisualizer):
    def __init__(self, singularity_graph):
        BlenderMeshVisualizer.__init__(self, singularity_graph.mesh)
        # singularity_graph
        # self.graph = HybridSingularityGraph(None, None)
        self.graph = singularity_graph

    # show base level : True
    # the input edge will be converted to background edge, curve is displayed.
    # show base level : False
    # the input edge does not include any curve information
    def show_singular_edge(self, singularity_edge_id,
                           radii=0.02, collection_name="Singular Edges",
                           rgba=None, show_ground_mesh_edge=True):
        if rgba is None:
            rgba = (1, 0, 0, 1)
        edge = self.graph.singular_Es[singularity_edge_id]
        if show_ground_mesh_edge:
            # edge to mesh eids
            self.show_edges(edge.es_link,
                            edge_names=[f'singular edge {singularity_edge_id}'] * len(edge.es_link),
                            radii=radii,
                            rgba=rgba,
                            show_ground_mesh_edge=show_ground_mesh_edge,
                            collection_name=collection_name)
            return
        vertex_list = self.graph.get_vertex_list()
        edge_vids = edge.start_end_vids
        if -1 in edge_vids:
            print(f'a circle singular edge {singularity_edge_id} is not displayed.')
            return
        make_edges(vertex_list, [edge_vids], r=radii,
                   name=collection_name, rgba=rgba, edge_names=[f'singular edge {singularity_edge_id}'])

    def show_singular_edges(self, eids,
                            radii=0.02, collection_name="Singular Edges",
                            rgba=None, show_ground_mesh_edge=True):
        generate_color_flag = False
        if rgba is None:
            generate_color_flag = True
        # different singular edge need assigned a new color
        for eid in eids:
            if generate_color_flag:
                pass
            self.show_singular_edge(eid,
                                    radii=radii,
                                    rgba=rgba,
                                    collection_name=collection_name,
                                    show_ground_mesh_edge=show_ground_mesh_edge)






