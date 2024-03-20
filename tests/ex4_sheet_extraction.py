from hybrid_meshio.hybrid_file_loader import MeshLoader
from os import path
from structure_analysis.sheet_extraction import SheetExtraction
from hybrid_base_complex.hybrid_singularity_structure import HybridSingularityGraph
from base_complex.base_complex import BaseComplex
from base_complex.bc_to_mesh_converter import BaseComplexToMeshConverter


def main_check_mesh():
    file_path = "../test_data/"
    file_name = "cube_twist-comp.obj.HYBRID"
    hybrid_mesh_loader = MeshLoader(file_path=path.join(file_path, file_name),
                                    load_as_hybrid=True)
    hybrid_mesh_loader.load()

    hsg = HybridSingularityGraph(hybrid_mesh_loader.mesh)
    bc = BaseComplex(hsg)
    bc.build()

    bcmc = BaseComplexToMeshConverter(bc)
    bcmc.covert()

    print(len(hybrid_mesh_loader.mesh.Cells))
    print("---------------------------")
    print(len(bcmc.new_mesh.Cells))

    se = SheetExtraction(bcmc.new_mesh)

    gs = se.mesh_cover_max_sheet_set()
    print(gs)


    return hybrid_mesh_loader.mesh


main_check_mesh()
