from pathlib import Path

from hybrid_meshio import hybrid_file_loader
from hybrid_base_complex import hybrid_singularity_structure
from hybrid_base_complex.hybrid_base_complex import BaseComplex
from base_complex.bc_to_mesh_converter import BaseComplexToMeshConverter
from structure_analysis.sheet_extraction import SheetExtraction
from structure_analysis.valence_singularity_graph_3d import ValenceSingularityGraph3D


def display_basic_mesh_property():
    file_dir = Path(__file__).resolve().parent.parent
    file_path = Path(file_dir).joinpath("test_data/")
    file_name = "cube_twist-comp.obj.HYBRID"

    hybrid_mesh_loader = hybrid_file_loader.MeshLoader(file_path=Path(file_path, file_name),
                                                       load_as_hybrid=True)

    hybrid_mesh_loader.load()

    hybrid_sg = hybrid_singularity_structure.HybridSingularityGraph(hybrid_mesh_loader.mesh)
    hybrid_sg.build()

    # base complex can accept hex dominant mesh by update the singular edges
    bc = BaseComplex(hybrid_sg, iteration_num=0)
    bc.build()

    bcmc = BaseComplexToMeshConverter(bc)
    bcmc.covert()

    se = SheetExtraction(bcmc.new_mesh)

    for sid, sheet in enumerate(se.sheet_objs):
        ssg_3d = ValenceSingularityGraph3D(se.mesh,
                                           parallel_eids=sheet.parallel_eids,
                                           region_cids=sheet.cells)
        ssg_3d.extract_valence_singular_graph(clean_graph=True)


display_basic_mesh_property()
