from hybrid_meshio.hybrid_file_loader import MeshLoader
from os import path


def main_check_mesh():
    file_path = "../test_data/"
    file_name = "twistcube_s.mesh"
    hybrid_mesh_loader = MeshLoader(file_path=path.join(file_path, file_name),
                                    load_as_hybrid=False)
    hybrid_mesh_loader.load()
    hybrid_mesh_loader.mesh.show_statistic()
    return hybrid_mesh_loader.mesh

main_check_mesh()
