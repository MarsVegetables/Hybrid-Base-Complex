from pathlib import Path

from hybrid_meshio import hybrid_file_loader
from hybrid_base_complex import hybrid_singularity_structure
from hybrid_base_complex.hybrid_base_complex import BaseComplex
from base_complex.bc_to_mesh_converter import BaseComplexToMeshConverter
from structure_analysis.sheet_extraction import SheetExtraction
from structure_analysis.sheet_max_matching import MaximumMatchingFinder
from structure_analysis.region_face_tracer import RegionFaceTracer
from structure_analysis.regional_singularity_structure_2d import SingularityStructure2D


def display_basic_mesh_property():
    file_dir = Path(__file__).resolve().parent.parent
    file_path = Path(file_dir).joinpath("test_data/")
    file_name = "twistcube_s.mesh"

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
    # ei.find_edges_not_in_sheets()
    imperfect_sids = [s.sheet_id for s in se.sheet_objs if not s.is_perfect]
    for i in imperfect_sids:
        intersecting_cids = se.sheet_objs[i].self_intersecting_cids
        if intersecting_cids:
            mmf = MaximumMatchingFinder(se.mesh, se.sheet_objs[i].parallel_eids, se.sheet_objs[i].cells)
            ams, amc = mmf.decompose_self_intersecting_sheet_layer()
            for j, ms in enumerate(ams):
                mc = amc[j]

                # face patches in the region
                rft = RegionFaceTracer(se.mesh, ms, mc)

                for cc, fp in enumerate(rft.face_patches):
                    ss_2d = SingularityStructure2D(se.mesh, mc, fp, ms)
                    ss_2d.build_singularity_graph()

                    # for jj, se_2d in enumerate(ss_2d.singular_es):
                    #     print(jj)


display_basic_mesh_property()
