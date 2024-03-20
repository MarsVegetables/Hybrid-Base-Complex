import sys
import os
from pathlib import Path

import bpy

# to load model in same dir
file_dir = Path(bpy.data.filepath).resolve().parent.parent
file_dir = f'{file_dir}{os.sep}'
if file_dir not in sys.path:
    sys.path.append(file_dir)

from hybrid_meshio import hybrid_file_loader
from blender_viz import mesh_visualizer
from hybrid_base_complex import hybrid_singularity_structure
from hybrid_base_complex.hybrid_base_complex import BaseComplex
from base_complex.bc_to_mesh_converter import BaseComplexToMeshConverter
from structure_analysis.sheet_extraction import SheetExtraction
from structure_analysis.sheet_max_matching import MaximumMatchingFinder
from structure_analysis.valence_singularity_graph_3d import ValenceSingularityGraph3D
from blender_viz.blender_viz_basic_op import color_id_to_rgba


def display_basic_mesh_property():
    file_path = None
    if not file_path:
        file_path = file_dir + "test_data/"
        print("file path is : ")
        print(file_path)
    file_name = None
    if not file_name:
        file_name = "cube_twist-comp.obj.HYBRID"

    hybrid_mesh_loader = hybrid_file_loader.MeshLoader(file_path=Path(file_path, file_name),
                                                       load_as_hybrid=True)

    hybrid_mesh_loader.load()

    bmv = mesh_visualizer.BlenderMeshVisualizer(hybrid_mesh_loader.mesh)
    bmv.show_wireframe()
    bmv.show_transparent_boundary(rgba=(0.2, 0.2, 0.2, 0.1))

    hybrid_sg = hybrid_singularity_structure.HybridSingularityGraph(hybrid_mesh_loader.mesh)
    hybrid_sg.build()

    # base complex can accept hex dominant mesh by update the singular edges
    bc = BaseComplex(hybrid_sg, iteration_num=0)
    bc.build()

    bcmc = BaseComplexToMeshConverter(bc)
    bcmc.covert()

    # visualize it
    bcmv = mesh_visualizer.BlenderMeshVisualizer(bcmc.new_mesh)

    se = SheetExtraction(bcmc.new_mesh)

    imperfect_sids = [s.sheet_id for s in se.sheet_objs if not s.is_perfect]
    for i in imperfect_sids:
        intersecting_cids = se.sheet_objs[i].self_intersecting_cids
        if intersecting_cids:
            bcmv.show_cells(cids=intersecting_cids, collection_name=f'interscting {i} cell')
            bcmv.show_cells(cids=se.sheet_objs[i].cells, collection_name=f'sheet {i} cell')
            mmf = MaximumMatchingFinder(se.mesh, se.sheet_objs[i].parallel_eids, se.sheet_objs[i].cells)
            ams, amc = mmf.decompose_self_intersecting_sheet_layer()
            for j, ms in enumerate(ams):
                mc = amc[j]
                bcmv.show_edges(eids=ms, collection_name=f"layer set {j} edges")
                bcmv.show_cells(cids=mc, collection_name=f"layer set {j} cells", alpha=0.1)

                ssg_3d = ValenceSingularityGraph3D(se.mesh, parallel_eids=ms, region_cids=mc)
                ssg_3d.extract_valence_singular_graph(clean_graph=True)
                # bcmv.show_edges(eids=ssg_3d.seed_singular_edges, collection_name=f"layer set {j} seed edges")
                for ii, vse in enumerate(ssg_3d.singular_es):
                    cell_rgba = color_id_to_rgba(color_id=vse.color,
                                                 total_color=ssg_3d.edge_color_num, alpha=1)
                    bcmv.show_edges(vse.es_link, collection_name=f"new se 11 set {j}", rgba=cell_rgba, radii=bcmv.default_r)


display_basic_mesh_property()
