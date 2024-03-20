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
from structure_analysis.regional_singularity_structure_2d import SingularityStructure2D
from blender_viz.blender_viz_basic_op import make_poly_line
from blender_viz.blender_viz_basic_op import color_id_to_rgba

from structure_analysis.region_face_tracer import RegionFaceTracer


def display_basic_mesh_property():
    file_path = None
    if not file_path:
        file_path = file_dir + "test_data/"
        print("file path is : ")
        print(file_path)
    file_name = None
    if not file_name:
        file_name = "twistcube_s.mesh"
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
    # ei.find_edges_not_in_sheets()
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
                bcmv.show_cells(cids=mc, collection_name=f"layer set {j} cells")

                # face patches in the region
                rft = RegionFaceTracer(se.mesh, ms, mc)
                for k, patch in enumerate(rft.face_patches):
                    bcmv.show_face(fids=patch, collection_name=f"s_{i}_m_{j}_patch_{k}")

                for cc, fp in enumerate(rft.face_patches):
                    ss_2d = SingularityStructure2D(se.mesh, mc, fp, ms)
                    ss_2d.build_singularity_graph()

                    bcmv.show_face(fids=ss_2d.face_scope, collection_name=f"s_{i}_m_{j}_face_scope_{cc}")

                    bcmv.show_vertices(vids=ss_2d.singular_vs, collection_name=f's_{i}_m_{j}_p_{cc} singular vert')
                    for jj, se_2d in enumerate(ss_2d.singular_es):
                        cell_rgba = color_id_to_rgba(color_id=se_2d.color,
                                                     total_color=ss_2d.color_num, alpha=1)
                        v = [ss_2d.mesh.Vertices[vid].xyz() for vid in se_2d.vs_link]
                        make_poly_line(v_list=v, obj_name=f"s_{i}_m_{j}_se_{jj}", collection_name=f"s_{i}_m_{j}",
                                       circle=se_2d.is_circle, rgba=cell_rgba, radius=bcmv.default_r)


display_basic_mesh_property()
