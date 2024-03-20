from hybrid_meshio.hybrid_file_loader import MeshLoader
from os import path
from structure_analysis.sheet_extraction import SheetExtraction
from hybrid_base_complex.hybrid_singularity_structure import HybridSingularityGraph
from base_complex.base_complex import BaseComplex
from base_complex.bc_to_mesh_converter import BaseComplexToMeshConverter
from structure_analysis.regional_singularity_structure_2d import SingularityStructure2D


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

    se = SheetExtraction(bcmc.new_mesh)

    for s in se.sheet_objs:
        if s.sheet_id != 0:
            continue
        for fids in s.separate_wall_fids:
            ss_2d = SingularityStructure2D(se.mesh, fids)
            ss_2d.build_singularity_graph()
            print(ss_2d.singular_vs)
            print(len(ss_2d.singular_es))
            # for es in ss_2d.singular_es[14]:
            ss_es = ss_2d.singular_es[14]
            print(f"es : {ss_es.es_link}")
            print(f"vs : {ss_es.vs_link}")
        break

    return hybrid_mesh_loader.mesh


main_check_mesh()
