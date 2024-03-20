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
from blender_viz import mesh_visualizer, singularity_graph_visualizer
from hybrid_base_complex import hybrid_singularity_structure


def display_basic_mesh_property():
    file_path = "../test_data/"
    file_name = "twistcube_s.mesh"
    hybrid_mesh_loader = hybrid_file_loader.MeshLoader(file_path=Path(file_path, file_name),
                                                       load_as_hybrid=True)

    hybrid_mesh_loader.load()

    bmv = mesh_visualizer.BlenderMeshVisualizer(hybrid_mesh_loader.mesh)
    bmv.show_wireframe()

    hybrid_sg = hybrid_singularity_structure.HybridSingularityGraph(hybrid_mesh_loader.mesh)
    hybrid_sg.build()

    sv = singularity_graph_visualizer.SingularityGraphVisualizer(hybrid_sg)
    # show singular edge
    # sv.show_singular_edge(0)
    # show multiple edge
    sv.show_singular_edges([i for i in range(len(hybrid_sg.singular_Es))], show_ground_mesh_edge=False)


display_basic_mesh_property()
