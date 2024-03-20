import sys
import os
from pathlib import Path

if "bpy" in locals():
    # not working
    # this next part forces a reload in case you edit the source after you first start the blender session
    import importlib

    importlib.reload(hybrid_meshio)
    importlib.reload(blender_viz)
else:
    import bpy

    # to load model in same dir
    file_dir = Path(bpy.data.filepath).resolve().parent.parent
    file_dir = f'{file_dir}{os.sep}'
    
    if file_dir not in sys.path:
        sys.path.append(file_dir)

    from hybrid_meshio import hybrid_file_loader
    from blender_viz import mesh_visualizer


def display_basic_mesh_property():
    file_path = input("please give a file path : ")
    if not file_path:
        file_path = file_dir + "test_data/"
        print("file path is : ")
        print(file_path)
    file_name = input("please give a file name : ")
    if not file_name:
        file_name = "cube_twist-comp.obj.HYBRID"

    hybrid_mesh_loader = hybrid_file_loader.MeshLoader(file_path=Path(file_path, file_name),
                                                       load_as_hybrid=True)
    
    hybrid_mesh_loader.load()

    bmv = mesh_visualizer.BlenderMeshVisualizer(hybrid_mesh_loader.mesh)
    bmv.show_wireframe()
    bmv.show_cells(cids=hybrid_mesh_loader.mesh.Cells)


display_basic_mesh_property()
