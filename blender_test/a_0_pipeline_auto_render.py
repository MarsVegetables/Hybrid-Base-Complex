import sys
import os
from pathlib import Path

import bpy

# to load model in same dir
file_dir = Path(bpy.data.filepath).resolve().parent.parent
file_dir = f'{file_dir}{os.sep}'

if file_dir not in sys.path:
    sys.path.append(file_dir)

from blender_viz.blender_auto_capture import CameraCapture
from blender_viz.blender_viz_basic_op import remove_collection, collection_render_control

from hybrid_meshio import hybrid_file_loader
from blender_viz import mesh_visualizer

from hybrid_base_complex.hybrid_singularity_structure import HybridSingularityGraph
from base_complex.base_complex import BaseComplex
from base_complex.bc_to_mesh_converter import BaseComplexToMeshConverter

from structure_analysis.valence_singularity_graph_3d import ValenceSingularityGraph3D
from blender_viz.blender_viz_basic_op import color_id_to_rgba

from structure_analysis.sheet_extraction import SheetExtraction
from structure_analysis.sheet_max_matching import MaximumMatchingFinder


def render_wireframe_and_face(capture_obj, target_mesh, img_name):
    face_collection_name = "mesh faces"
    bmv = mesh_visualizer.BlenderMeshVisualizer(target_mesh)
    wireframe_collection_name = "mesh wireframe"
    wireframe_obj_name = bmv.show_wireframe(collection_name=wireframe_collection_name)
    bmv.show_face(fids=target_mesh.Faces, face_name="all_faces", collection_name=face_collection_name)

    # save
    capture_obj.update_camera_and_light(target_name=wireframe_obj_name)
    capture_obj.render_view(img_name=img_name, target_info=target_mesh.__str__())
    # remove_collection(wireframe_collection_name)
    collection_render_control(collection_name=wireframe_collection_name, hide_render=True)
    remove_collection(face_collection_name)
    return wireframe_obj_name, wireframe_collection_name, face_collection_name


def show_singularity_edges(singular_graph, visualizer, collection_name="singularity edges",
                           show_redundant=False, show_color=False, adjust_alpha=False):
    singular_objs = singular_graph.singular_es
    edge_rgba = (1, 0, 0, 1)
    if show_redundant:
        singular_objs = singular_graph.redundant_singular_es
        edge_rgba = (0, 0.8, 0.2, 1)
    for i, se in enumerate(singular_objs):
        if show_color:
            edge_rgba = color_id_to_rgba(color_id=se.color,
                                         total_color=singular_graph.edge_color_num, alpha=1)
        edge_end_vids = set()
        if adjust_alpha and not se.is_circle:
            for vid in se.start_end_vids:
                # redundant vert does not impact valence
                if vid in singular_graph.redundant_singular_vs:
                    continue
                edge_end_vids.add(vid)
        visualizer.show_edges(se.es_link, collection_name=collection_name,
                              rgba=edge_rgba, colored_vids=edge_end_vids,
                              adjust_alpha=adjust_alpha)


# def show_region_faces(visualizer):
#     # face patches in the region
#     rft = RegionFaceTracer(se.mesh, ms, mc)
#     for k, patch in enumerate(rft.face_patches):
#         visualizer.show_face(fids=patch, collection_name=f"s_{i}_m_{j}_patch_{k}")

def valence_based_singularity_graph(mesh, cids, sheet_parallel_eids,
                                    capture_obj, img_section_name,
                                    remove_parallel_edge_in_sg=True):
    bcmv = mesh_visualizer.BlenderMeshVisualizer(mesh)

    # valence-based singularity graph
    ssg_3d = ValenceSingularityGraph3D(mesh, parallel_eids=sheet_parallel_eids,
                                       region_cids=cids)

    edge_collection_name = "singularity edges"

    # init edges
    ssg_3d.init_structure_singularity_graph()
    bcmv.show_edges(eids=ssg_3d.seed_singular_edges, collection_name=edge_collection_name)
    img_name = f"{img_section_name}_t_0_valence_seed_se"
    capture_obj.render_view(img_name=img_name, target_info=mesh.__str__())
    remove_collection(edge_collection_name)

    # all init se
    show_singularity_edges(ssg_3d, bcmv, collection_name=edge_collection_name)
    img_name = f"{img_section_name}_t_1_valence_all_init_se"
    capture_obj.render_view(img_name=img_name, target_info=mesh.__str__())
    remove_collection(edge_collection_name)

    # complete
    ssg_3d.complete_valence_singularity_graph()
    show_singularity_edges(ssg_3d, bcmv, collection_name=edge_collection_name)
    img_name = f"{img_section_name}_t_2_complete_valence_se"
    capture_obj.render_view(img_name=img_name, target_info=mesh.__str__())
    remove_collection(edge_collection_name)

    # complete valence singularity edges
    ssg_3d.clean_valence_singularity_graph(remove_parallel_edge_in_sg)
    ssg_3d.final_singular_elements()

    # show singular edges
    show_singularity_edges(ssg_3d, bcmv, show_color=True,
                           collection_name=edge_collection_name, adjust_alpha=True)
    img_name = f"{img_section_name}_t_3_cleaned_valence_se"
    capture_obj.render_view(img_name=img_name, target_info=mesh.__str__())
    remove_collection(edge_collection_name)

    # show redundant singular edges
    show_singularity_edges(ssg_3d, bcmv, show_redundant=True,
                           collection_name=edge_collection_name)
    img_name = f"{img_section_name}_t_4_redundant_valence_se"
    capture_obj.render_view(img_name=img_name, target_info=mesh.__str__())
    remove_collection(edge_collection_name)


def mesh_structure_visualization_pipeline(file_folder, file_name):
    hybrid_mesh_loader = hybrid_file_loader.MeshLoader(file_path=Path(file_folder, file_name),
                                                       load_as_hybrid=True)

    if not hybrid_mesh_loader.load():
        return

    # create camera capture object
    # specific the file path and file name
    # the capture obj will create a render folder under the file path
    # and a fold with file_name will also create in render folder
    cc = CameraCapture(file_folder, file_name, show_mesh_info=False)

    # target_edge_length = cur_mesh.get_average_edge_length()
    print("step_0_mesh_surface")
    _, wireframe_collection_name, _ = render_wireframe_and_face(cc, hybrid_mesh_loader.mesh,
                                                                "step_0_mesh_surface")

    # mesh visualization
    bmv = mesh_visualizer.BlenderMeshVisualizer(hybrid_mesh_loader.mesh)
    # wireframe_collection_name = "mesh wireframe"
    # wireframe_obj_name = bmv.show_wireframe(collection_name=wireframe_collection_name)

    # show transparent bnd
    print("show transparent bnd")
    transparent_bnd_collection_name = "trans bnd"
    bmv.show_transparent_boundary(collection_name=transparent_bnd_collection_name)
    # collection_render_control(collection_name=transparent_bnd_collection_name, hide_render=True)
    # remove_collection(transparent_bnd_collection_name)

    # hybrid singularity graph
    print("step_1_hybrid_singularity_graph")
    hsg = HybridSingularityGraph(hybrid_mesh_loader.mesh)
    init_hybrid_sg_collection_name = "init hybrid sg"
    for se in hsg.singular_Es:
        bmv.show_edges(eids=se.es_link, collection_name=init_hybrid_sg_collection_name)
    cc.render_view(img_name=f"step_1_hybrid_singularity_graph", target_info=hybrid_mesh_loader.mesh.__str__())
    remove_collection(init_hybrid_sg_collection_name)

    # base complex
    bc = BaseComplex(hsg)
    bc.build()

    bcmc = BaseComplexToMeshConverter(bc)
    bcmc.covert()

    bcmv = mesh_visualizer.BlenderMeshVisualizer(bcmc.new_mesh)

    print("step_2_complete_hybrid_singularity_graph")
    # complete hybrid singularity graph
    complete_hybrid_sg_collection_name = "complete hybrid sg"
    bcmv.show_edges(eids=bcmc.new_mesh.Edges, collection_name=complete_hybrid_sg_collection_name)
    # render
    cc.render_view(img_name=f"step_2_complete_hybrid_singularity_graph", target_info=hybrid_mesh_loader.mesh.__str__())
    remove_collection(complete_hybrid_sg_collection_name)

    # hybrid base complex
    print("step_3_hybrid_base_complex")
    # bcmv.show_cells(cids=bcmc.new_mesh.get_non_hexa_cell_id_list(), collection_name="non hex")
    hybrid_bc_collection_name = "hybrid base complex"
    bcmv.show_cells(cids=bcmc.new_mesh.Cells, collection_name=hybrid_bc_collection_name)
    # render
    collection_render_control(collection_name=transparent_bnd_collection_name, hide_render=True)
    collection_render_control(collection_name=wireframe_collection_name, hide_render=False)
    cc.render_view(img_name=f"step_3_hybrid_base_complex", target_info=hybrid_mesh_loader.mesh.__str__())
    remove_collection(hybrid_bc_collection_name)
    # update surface and wireframe visibility
    collection_render_control(collection_name=wireframe_collection_name, hide_render=True)
    collection_render_control(collection_name=transparent_bnd_collection_name, hide_render=False)

    # mesh valence-based singularity graph
    print("step_4_mesh")
    valence_based_singularity_graph(mesh=bcmc.new_mesh,
                                    cids=set(),
                                    sheet_parallel_eids=set(),
                                    capture_obj=cc,
                                    img_section_name="step_4_mesh")
    print("step_5_largest_sheet")
    parallel_edge_color = (0.3, 0.3, 0.8, 1)
    sheet_edge_color = (0.1, 0.1, 0.1, 1)
    # sheets
    sheet_extractor = SheetExtraction(bcmc.new_mesh)
    # largest sheets
    sheet_ids = sheet_extractor.get_largest_sheet_ids()
    # valence-based singularity graph
    for i in sheet_ids:
        sheet = sheet_extractor.sheet_objs[i]
        # render sheet
        sheet_cell_collection_name = "sheet elements"
        sheet_eids = set()
        for cid in sheet.cells:
            cell = bcmc.new_mesh.Cells[cid]
            sheet_eids.update(cell.Eids)
        sheet_eids = sheet_eids.difference(sheet.parallel_eids)
        bcmv.show_edges(eids=sheet_eids, collection_name=sheet_cell_collection_name,
                        radii=bcmv.default_r * 0.2, rgba=sheet_edge_color)
        bcmv.show_edges(eids=sheet.parallel_eids, collection_name=sheet_cell_collection_name,
                        radii=bcmv.default_r * 0.5, rgba=parallel_edge_color)
        bcmv.show_cells(cids=sheet.cells, collection_name=sheet_cell_collection_name)
        cc.render_view(img_name=f"step_5_largest_sheet_{i}", target_info=hybrid_mesh_loader.mesh.__str__())
        remove_collection(sheet_cell_collection_name)
        # valence sg
        valence_based_singularity_graph(mesh=bcmc.new_mesh,
                                        cids=sheet.cells,
                                        sheet_parallel_eids=sheet.parallel_eids,
                                        capture_obj=cc,
                                        img_section_name=f"step_5_largest_sheet_{i}",
                                        remove_parallel_edge_in_sg=False)

    # imperfect sheets
    print("step_6_imperfect_sheet")
    # imperfect_sids = sheet_extractor.get_imperfect_sheet_ids()
    ranked_imperfect_sids = sheet_extractor.get_size_ranked_imperfect_sheet_ids()
    # valence-based singularity graph
    for imperfect_sid in ranked_imperfect_sids[:5]:
        sheet = sheet_extractor.sheet_objs[imperfect_sid]
        # render sheet
        sheet_cell_collection_name = "sheet elements"
        number_of_cells = len(sheet.cells)
        sheet_eids = set()
        for cid in sheet.cells:
            cell = bcmc.new_mesh.Cells[cid]
            sheet_eids.update(cell.Eids)
        sheet_eids = sheet_eids.difference(sheet.parallel_eids)
        bcmv.show_edges(eids=sheet_eids, collection_name=sheet_cell_collection_name,
                        radii=bcmv.default_r*0.2, rgba=sheet_edge_color)
        bcmv.show_edges(eids=sheet.parallel_eids, collection_name=sheet_cell_collection_name,
                        radii=bcmv.default_r * 0.5, rgba=parallel_edge_color)
        bcmv.show_cells(cids=sheet.cells, collection_name=sheet_cell_collection_name)
        cc.render_view(img_name=f"step_6_sheet_{imperfect_sid}_imperfect_cell_{number_of_cells}",
                       target_info=hybrid_mesh_loader.mesh.__str__())
        remove_collection(sheet_cell_collection_name)
        # render sheet valence sg
        valence_based_singularity_graph(mesh=bcmc.new_mesh,
                                        cids=sheet.cells,
                                        sheet_parallel_eids=sheet.parallel_eids,
                                        capture_obj=cc,
                                        img_section_name=f"step_6_sheet_{imperfect_sid}_imperfect_cell_{number_of_cells}",
                                        remove_parallel_edge_in_sg=False)

        # step 6 self-intersecting sheet decomposition and visualization
        intersecting_cids = sheet.self_intersecting_cids

        if not intersecting_cids:
            continue

        mmf = MaximumMatchingFinder(sheet_extractor.mesh,
                                    sheet.parallel_eids,
                                    sheet.cells)
        ams, amc = mmf.decompose_self_intersecting_sheet_layer()
        for j, ms in enumerate(ams):
            mc = amc[j]
            # render sheet
            sheet_cell_collection_name = "sheet elements"
            img_section_name = f"step_6_sheet_{imperfect_sid}_self_intersecting_subset_{j}"
            sheet_eids = set()
            for cid in mc:
                cell = bcmc.new_mesh.Cells[cid]
                sheet_eids.update(cell.Eids)
            sheet_eids = sheet_eids.difference(ms)
            bcmv.show_edges(eids=sheet_eids, collection_name=sheet_cell_collection_name,
                            radii=bcmv.default_r * 0.2, rgba=sheet_edge_color)
            bcmv.show_edges(eids=ms, collection_name=sheet_cell_collection_name,
                            radii=bcmv.default_r * 0.5, rgba=parallel_edge_color)
            bcmv.show_cells(cids=mc, collection_name=sheet_cell_collection_name)
            cc.render_view(img_name=img_section_name,
                           target_info=hybrid_mesh_loader.mesh.__str__())
            remove_collection(sheet_cell_collection_name)
            # subset valence-based sg
            valence_based_singularity_graph(mesh=bcmc.new_mesh,
                                            cids=mc,
                                            sheet_parallel_eids=ms,
                                            capture_obj=cc,
                                            img_section_name=img_section_name)

    # clean
    print("step_7_cleaning")
    cc.clean_old_collection(reset_img_idx=True)


def find_all_hex_dom_files(dirPath):
    hybrid_files = []
    for root, dirs, d_files in os.walk(dirPath):
        for file in d_files:
            if file.endswith(".HYBRID") or file.endswith(".mesh") or file.endswith(".hedra"):
                # hybrid_files.append(os.path.join(root, file))
                hybrid_files.append(file)
    return hybrid_files


if __name__ == '__main__':
    file_path = Path(file_dir + "test_data")
    files = find_all_hex_dom_files(file_path)
    for file_n in files:
        mesh_structure_visualization_pipeline(file_path, file_n)
