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


def show_mesh_sheets():
    file_path = input("please give a file path : ")
    if not file_path:
        file_path = file_dir + "test_data/"
        print("file path is : ")
        print(file_path)
    file_name = input("please give a file name : ")
    if not file_name:
        file_name = "cylinder_n.mesh"

    hybrid_mesh_loader = MeshLoader(file_path=Path(file_path, file_name),
                                    load_as_hybrid=True)
    hybrid_mesh_loader.load()
    print(hybrid_mesh_loader.mesh.mesh_level)

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
    for sheet in se.sheet_objs:
        bcmv.show_cells(cids=sheet.cells, face_name=f'cell', collection_name=f'sheet-{sheet.sheet_id}')
        bcmv.show_cells(cids=sheet.non_hexa_cids, collection_name=f'sheet {sheet.sheet_id} non hex cell')

    # return hybrid_mesh_loader.mesh


show_mesh_sheets()
