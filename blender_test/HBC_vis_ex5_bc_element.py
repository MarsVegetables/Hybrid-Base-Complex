import sys
import os
from pathlib import Path

import bpy

# to load model in same dir
file_dir = Path(bpy.data.filepath).resolve().parent.parent
file_dir = f'{file_dir}{os.sep}'
if file_dir not in sys.path:
    sys.path.append(file_dir)

from blender_viz import mesh_visualizer

from hybrid_meshio.hybrid_file_loader import MeshLoader
from structure_analysis.sheet_extraction import SheetExtraction
from hybrid_base_complex.hybrid_singularity_structure import HybridSingularityGraph
from base_complex.base_complex import BaseComplex
from base_complex.bc_to_mesh_converter import BaseComplexToMeshConverter
from structure_analysis.valence_singularity_graph_3d import ValenceSingularityGraph3D
from blender_viz.blender_viz_basic_op import color_id_to_rgba


def show_mesh_sheets():
    file_path = input("please give a file path : ")
    if not file_path:
        file_path = file_dir + "test_data/"
        print("file path is : ")
        print(file_path)
    file_name = input("please give a file name : ")
    if not file_name:
        file_name = "test_2.HYBRID"
    hybrid_mesh_loader = MeshLoader(file_path=Path(file_path, file_name),
                                    load_as_hybrid=True)
    hybrid_mesh_loader.load()

    hsg = HybridSingularityGraph(hybrid_mesh_loader.mesh)
    bc = BaseComplex(hsg)
    bc.build()

    bcmc = BaseComplexToMeshConverter(bc)
    bcmc.covert()

    se = SheetExtraction(bcmc.new_mesh)

    bmv = mesh_visualizer.BlenderMeshVisualizer(hybrid_mesh_loader.mesh)
    bmv.show_wireframe()

    # visualize it
    bcmv = mesh_visualizer.BlenderMeshVisualizer(bcmc.new_mesh)
    bcmv.show_cells(cids=bcmc.new_mesh.get_non_hexa_cell_id_list(), collection_name="non hex")

    # bcmv.show_edges(eids=bcmc.new_mesh.Edges, collection_name="all HBC edges")

    ssg_3d = ValenceSingularityGraph3D(se.mesh)
    ssg_3d.init_structure_singularity_graph()
    bcmv.show_edges(eids=ssg_3d.seed_singular_edges, collection_name=f"seed edges")
    bcmv.show_vertices(vids=ssg_3d.singular_vs, collection_name="sv")
    ssg_3d.complete_seed_edges()
    bcmv.show_edges(eids=ssg_3d.seed_singular_edges, collection_name=f"seed edges after")
    ssg_3d.complete_singular_vert()
    bcmv.show_vertices(vids=ssg_3d.singular_vs, collection_name="sv after")
    ssg_3d.complete_valence_singularity_graph()
    print("render edges")
#    for i, se in enumerate(ssg_3d.singular_es):
#        bcmv.show_edges(se.es_link, collection_name=f"se {i}")

    ssg_3d.clean_valence_singularity_graph()
    for i, se in enumerate(ssg_3d.singular_es):
        cell_rgba = color_id_to_rgba(color_id=se.color,
                                     total_color=ssg_3d.edge_color_num, alpha=1)
        bcmv.show_edges(se.es_link, collection_name=f"new se", rgba=cell_rgba, radii=bcmv.default_r)
    for i, se in enumerate(ssg_3d.redundant_singular_es):
        bcmv.show_edges(se.es_link, collection_name=f"new re se", rgba=(0, 0, 1, 1),
                        radii=bcmv.default_r)
    ssg_3d.final_singular_elements()
    for i, se in enumerate(ssg_3d.singular_es):
        cell_rgba = color_id_to_rgba(color_id=se.color,
                                     total_color=ssg_3d.edge_color_num, alpha=1)
        bcmv.show_edges(se.es_link, collection_name=f"new se 11", rgba=cell_rgba, radii=bcmv.default_r)
    for i, se in enumerate(ssg_3d.redundant_singular_es):
        bcmv.show_edges(se.es_link, collection_name=f"new re se 11", rgba=(0, 0, 1, 1),
                        radii=bcmv.default_r)

    # return hybrid_mesh_loader.mesh


show_mesh_sheets()
