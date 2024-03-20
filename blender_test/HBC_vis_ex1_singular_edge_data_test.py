from pathlib import Path
from hybrid_meshio.hybrid_file_loader import MeshLoader

from hybrid_base_complex.hybrid_singularity_structure import HybridSingularityGraph


def display_basic_mesh_property():
    file_path = "../test_data/"
    file_name = "twistcube_s.mesh"
    hybrid_mesh_loader = MeshLoader(file_path=Path(file_path, file_name),
                                    load_as_hybrid=True)

    hybrid_mesh_loader.load()

    hybrid_sg = HybridSingularityGraph(hybrid_mesh_loader.mesh)
    hybrid_sg.build()

    edge = hybrid_sg.singular_Es[21]

    (current_mesh,
     edge_vids,
     current_edge_names) = hybrid_sg.mesh.covert_to_background_edge_vids(edge.es_link,
                                                                         edge_names=[f'singular edge {0}'] * len(edge.es_link),
                                                                         show_ground_mesh_edge=False)
    print(current_mesh)
    print(edge_vids)
    print(current_edge_names)



display_basic_mesh_property()
