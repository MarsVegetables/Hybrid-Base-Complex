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


def display_basic_mesh_property():
    file_path = file_dir + "test_data/"
    file_name = "teapot.hedra"

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
    bcmv.show_cells(cids=bcmv.mesh.Cells)


display_basic_mesh_property()
