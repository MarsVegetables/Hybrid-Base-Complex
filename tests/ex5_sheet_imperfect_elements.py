from os import path
from pathlib import Path
from hybrid_meshio.hybrid_file_loader import MeshLoader
from hybrid_base_complex.hybrid_singularity_structure import HybridSingularityGraph
from base_complex.base_complex import BaseComplex
from base_complex.bc_to_mesh_converter import BaseComplexToMeshConverter
from structure_analysis.sheet_identification import SheetIdentification


def main_check_mesh():
    file_path = input("please give a file path : ")
    if not file_path:
        file_dir = Path(__file__).resolve().parent.parent
        file_path = file_dir.joinpath("test_data/")
        print("file path is : ")
        print(file_path)
    file_name = input("please give a file name : ")
    if not file_name:
        file_name = "hanger_output.HYBRID"
    hybrid_mesh_loader = MeshLoader(file_path=path.join(file_path, file_name),
                                    load_as_hybrid=True)
    hybrid_mesh_loader.load()
    hsg = HybridSingularityGraph(hybrid_mesh_loader.mesh)
    bc = BaseComplex(hsg)
    bc.build()

    bcmc = BaseComplexToMeshConverter(bc)
    bcmc.covert()

    se = SheetIdentification(bcmc.new_mesh)
    for s in se.sheet_objs:
        if s.is_perfect:
            continue
        c, e, v = se.get_sheet_imperfect_element_ids(s.sheet_id)
        print("----------------------------------------------------------------")
        print(c)
        print(e)
        print(v)


    return bcmc.new_mesh


main_check_mesh()
