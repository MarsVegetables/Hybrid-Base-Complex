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
from blender_viz.blender_viz_basic_op import color_id_to_rgba
from structure_analysis.sheet_identification import SheetIdentification


def display_basic_mesh_property():
    file_path = input("please give a file path : ")
    if not file_path:
        file_path = file_dir + "test_data/"
        print("file path is : ")
        print(file_path)
    file_name = input("please give a file name : ")
    if not file_name:
        file_name = "hanger_output.HYBRID"

    hybrid_mesh_loader = hybrid_file_loader.MeshLoader(file_path=Path(file_path, file_name),
                                                       load_as_hybrid=True)

    hybrid_mesh_loader.load()

    bmv = mesh_visualizer.BlenderMeshVisualizer(hybrid_mesh_loader.mesh)
    bmv.show_wireframe()

    hybrid_sg = hybrid_singularity_structure.HybridSingularityGraph(hybrid_mesh_loader.mesh)
    hybrid_sg.build()

    # base complex can accept hex dominant mesh by update the singular edges
    bc = BaseComplex(hybrid_sg, iteration_num=0)
    bc.build()

    bcmc = BaseComplexToMeshConverter(bc)
    bcmc.covert()
    
    # visualize it
    bcmv = mesh_visualizer.BlenderMeshVisualizer(bcmc.new_mesh)

    ei = SheetIdentification(bcmc.new_mesh)

    unmatched_vids = [v for v in ei.vert_weight if ei.vert_weight[v] > 0]
    bcmv.show_vertices(vids=unmatched_vids, collection_name=f'all vert', radii_weight=ei.vert_weight)

    imperfect_sids = [s.sheet_id for s in ei.sheet_objs if not s.is_perfect]
    for i in imperfect_sids:
        sheet = ei.sheet_objs[i]
        sheet_rgba = color_id_to_rgba(sheet.color_id, total_color=ei.color_num)

        bcmv.show_cells(cids=ei.sheet_objs[i].cells, rgba=sheet_rgba, collection_name=f'imperfect sheet {i} cell')
        bcmv.show_edges(eids=ei.sheet_objs[i].parallel_eids, collection_name=f'sheet {i} parallel edges')
        # element vis
        c, e, v = ei.get_sheet_imperfect_element_ids(i)
        if c:
            bcmv.show_cells(cids=c, rgba=sheet_rgba, collection_name=f'imperfect sheet {i} self intersecting cell')
        # if e:
        #     bcmv.show_edges(eids=e, collection_name=f'imperfect sheet {i} self parallel edges')
        if v:
            bcmv.show_vertices(vids=v, collection_name=f'imperfect sheet {i} ill vert')


display_basic_mesh_property()
