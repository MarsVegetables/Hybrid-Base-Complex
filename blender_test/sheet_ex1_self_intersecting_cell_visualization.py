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


def display_basic_mesh_property():
    file_path = None
    if not file_path:
        file_path = file_dir + "test_data/"
        print("file path is : ")
        print(file_path)
    file_name = input("please give a file name : ")
    if not file_name:
        file_name = "twistcube_s.mesh"
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

    se = SheetExtraction(bcmc.new_mesh)

    imperfect_sids = [s.sheet_id for s in se.sheet_objs if not s.is_perfect]
    for i in imperfect_sids:
        intersecting_cids = se.sheet_objs[i].self_intersecting_cids
        if intersecting_cids:
            bcmv.show_cells(cids=intersecting_cids, collection_name=f'interscting {i} cell')
            bcmv.show_cells(cids=se.sheet_objs[i].cells, collection_name=f'sheet {i} cell')




display_basic_mesh_property()
