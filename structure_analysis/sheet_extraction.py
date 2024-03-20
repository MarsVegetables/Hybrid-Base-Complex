from structure_analysis.sheet_elements import *
from mesh_structure.mesh_common import check_common_elements
from tqdm import tqdm
import random
import numpy as np


class SheetExtraction:
    def __init__(self, mesh):
        self.mesh = mesh
        # self.mesh = BasicMesh("", "")
        self.sheets_list = []
        self.sheet_objs = []
        self.max_neighbor_num = 0
        self.color_num = 0
        # temporary usage
        self.cell_visited_map = []
        # end
        self.cell_sheet_map = {cid: set() for cid in self.mesh.Cells}
        # edges that in sheets that contain at least one hex
        self.edge_sheet_map = {eid: set() for eid in self.mesh.Edges}
        self.build_sheets()
        self.build_sheet_attributes()
        self.build_sheet_neighbor_connectivity()

    # -------------- sheet building process -------------------
    # neighbors of a parallel edge are considered in a sheet
    def build_sheets(self):
        print("building sheets ...")
        # clean visit flag
        self.mesh.clean_visited()
        unvisted_edges = list(self.mesh.Edges)
        while unvisted_edges:
            eid = unvisted_edges.pop()
            self.build_sheet_by_edge(eid)
            if not self.mesh.Edges[eid].visited:
                unvisted_edges.append(eid)
        print("done ...")

    # split sheet
    def build_sheet_by_edge(self, eid):
        edge = self.mesh.Edges[eid]
        if edge.visited:
            return None
        sheet_cids, sheet_parallel_eids = self.trace_cell_and_parallel_edges(eid)
        if len(sheet_cids) == 0:
            return
        # split sheets
        (sheet_cids_groups,
         sheet_parallel_eids_groups) = self.split_sheet_by_cell_connection(sheet_cids, sheet_parallel_eids)
        for s, p in zip(sheet_cids_groups, sheet_parallel_eids_groups):
            if len(s) == 0:
                continue
            if self.is_sheet_exist(s, p):
                continue
            # build sheet object
            sheet = HybridSheet(sheet_id=len(self.sheet_objs))
            sheet.parallel_eids = p
            sheet.cells = s
            sheet.hexa_cids = set(s)
            self.sheet_objs.append(sheet)
            self.build_edge_sheet_map(sheet)
        return 0

    def is_sheet_exist(self, sheet_cids, sheet_parallel_eids):
        find_flag = False
        for sheet in self.sheet_objs:
            if sheet_cids != sheet.cells:
                continue
            if sheet.parallel_eids != sheet_parallel_eids:
                continue
            find_flag = True
        return find_flag

    def build_edge_sheet_map(self, sheet):
        if len(sheet.hexa_cids) == 0:
            return
        else:
            for eid in sheet.parallel_eids:
                self.edge_sheet_map[eid].add(sheet.sheet_id)

    # in sheet cells connect by edge or face
    def trace_cell_and_parallel_edges(self, eid):
        eid_stack = [eid]
        parallel_eids = set(eid_stack)
        sheet_cids = set()
        while eid_stack:
            cur_eid = eid_stack.pop()
            edge = self.mesh.Edges[cur_eid]
            if edge.visited:
                continue
            edge.visited = True
            # only contain hex
            cell_eids = set()
            for cid in edge.neighboring_Cids:
                cell = self.mesh.Cells[cid]
                if cell.is_hexa():
                    sheet_cids.add(cid)
                    cell_eids.update(cell.Eids)
            next_parallel_edges = check_common_elements(edge.parallels, cell_eids)
            eid_stack.extend(next_parallel_edges)
            parallel_eids.update(next_parallel_edges)
        return list(sheet_cids), list(parallel_eids)

    def split_sheet_by_cell_connection(self, sheet_cids, parallel_eids):
        # self.mesh.clean_visited()
        self.cell_visited_map = [False] * len(self.mesh.Cells)
        connected_cid_groups = []
        for cid in sheet_cids:
            if self.cell_visited_map[cid]:
                continue
            cids = self.search_connected_cid_in_sheet_scope(cid, sheet_cids)
            connected_cid_groups.append(cids)
        if len(connected_cid_groups) <= 1:
            return connected_cid_groups, [parallel_eids]
        # eids
        corresponding_eids = []
        for cids in connected_cid_groups:
            group_parallel_eids = set()
            for cid in cids:
                cell = self.mesh.Cells[cid]
                group_parallel_eids.update(cell.Eids)
            group_parallel_eids = check_common_elements(group_parallel_eids, parallel_eids)
            corresponding_eids.append(group_parallel_eids)
        return connected_cid_groups, corresponding_eids

    def search_connected_cid_in_sheet_scope(self, cid, sheet_cids):
        bfs_stack = [cid]
        connected_cids = set()
        while bfs_stack:
            cur_cid = bfs_stack.pop()
            if self.cell_visited_map[cur_cid]:
                continue
            self.cell_visited_map[cur_cid] = True
            cell = self.mesh.Cells[cur_cid]
            connected_cids.add(cur_cid)
            cc = check_common_elements(sheet_cids, cell.neighboring_Cids)
            bfs_stack.extend(cc)
        return connected_cids

    # -------------- sheet building process done -------------------

    # -------------- sheet attributes building process -------------------
    def build_sheet_attributes(self):
        print("building sheet attributes ...")
        self.create_matching_map()
        self.add_non_hex_to_sheet()
        self.identify_sheet_cells_type()
        # self.build_edge_sheet_map()
        self.find_sheet_boundary_faces()
        print("done ...")

    def add_non_hex_to_sheet(self):
        for sheet in self.sheet_objs:
            for eid in sheet.parallel_eids:
                edge = self.mesh.Edges[eid]
                for cid in edge.neighboring_Cids:
                    if self.mesh.Cells[cid].is_hexa():
                        continue
                    sheet.non_hexa_cids.add(cid)

    def identify_sheet_cells_type(self):
        for sheet in self.sheet_objs:
            for cid in sheet.cells:
                if self.mesh.Cells[cid].is_hexa():
                    sheet.hexa_cids.add(cid)
                else:
                    sheet.non_hexa_cids.add(cid)

    # def build_edge_sheet_map(self, sheet):
    #     for sheet in self.sheet_objs:
    #         if len(sheet.hexa_cids) == 0:
    #             continue
    #         else:
    #             for eid in sheet.parallel_eids:
    #                 self.edge_sheet_map[eid].add(sheet.sheet_id)

    # if a cell's all neighboring cell are in same sheet
    # the cell is inside the sheet
    # only consider hex cells,
    # non hex cell will make sheet boundary very complex
    def is_cell_inside_sheet(self, sheet_id, cid):
        ncids = self.mesh.Cells[cid].neighboring_Cids
        tmp = check_common_elements(ncids, self.sheet_objs[sheet_id].hexa_cids)
        if len(ncids) == len(tmp):
            return True
        return False

    def find_sheet_boundary_faces(self):
        for sheet in self.sheet_objs:
            for cid in sheet.cells:
                if self.is_cell_inside_sheet(sheet.sheet_id, cid):
                    continue
                self.find_boundary_faces_of_cell_in_sheet(cid, sheet.sheet_id)
            self.identify_wall_face_of_sheet(sheet.sheet_id)
            self.split_wall_fids(sheet.sheet_id)

    def find_boundary_faces_of_cell_in_sheet(self, cid, sheet_id):
        for fid in self.mesh.Cells[cid].Fids:
            face = self.mesh.Faces[fid]
            sheet = self.sheet_objs[sheet_id]
            tmp = check_common_elements(face.neighboring_Cids,
                                        sheet.cells)
            if len(tmp) == len(face.neighboring_Cids) and len(face.neighboring_Cids) != 1:
                continue
            sheet.boundary_fids.add(fid)
            # a face can contain many neighboring sheet
            # in three different direction
            face.neighboring_sheet_ids.add(sheet_id)

    def identify_wall_face_of_sheet(self, sheet_id):
        sheet = self.sheet_objs[sheet_id]
        boundary_fids = list(sheet.boundary_fids)
        for fid in boundary_fids:
            cur_face = self.mesh.Faces[fid]
            eids = cur_face.Eids
            tmp = check_common_elements(eids, sheet.parallel_eids)
            if len(tmp) > 0:
                continue
            sheet.wall_fids.append(fid)

    def filter_out_share_cell_neighbor_faces(self, fid, wall_fids):
        flat_neighbor_fids = []
        face = self.mesh.Faces[fid]
        neighbor_fids = face.neighboring_Fids
        for nfid in neighbor_fids:
            if nfid not in wall_fids:
                continue
            n_face = self.mesh.Faces[nfid]
            tmp = check_common_elements(face.neighboring_Cids, n_face.neighboring_Cids)
            # contain same neighbor cell means two faces are in the same cell
            if len(tmp) > 0:
                continue
            flat_neighbor_fids.append(nfid)
        return flat_neighbor_fids

    def trace_wall_patch(self, fid, wall_fids, visited_flags):
        wall_patch = set()
        fid_stack = [fid]
        while fid_stack:
            cur_fid = fid_stack.pop()
            if visited_flags[cur_fid]:
                continue
            visited_flags[cur_fid] = True
            wall_patch.add(cur_fid)
            next_fids = self.filter_out_share_cell_neighbor_faces(cur_fid, wall_fids)
            fid_stack.extend(next_fids)
        return wall_patch

    def split_wall_fids(self, sheet_id):
        sheet = self.sheet_objs[sheet_id]
        visited_flag = {fid: False for fid in sheet.wall_fids}
        for wfid in sheet.wall_fids:
            if visited_flag[wfid]:
                continue
            wall_patch = self.trace_wall_patch(wfid, sheet.wall_fids, visited_flag)
            sheet.separate_wall_fids.append(wall_patch)

    # -------------- sheet attributes building process done -------------------

    # -------------- sheet relationship building process -------------------
    def find_sheet_neighbor_sheet(self, sheet):
        wall_neighbor_sheet_ids = {}
        for i, wall_fids in enumerate(sheet.separate_wall_fids):
            neighbor_sid = set()
            for fid in wall_fids:
                cur_face = self.mesh.Faces[fid]
                for sid in cur_face.neighboring_sheet_ids:
                    if sid == sheet.sheet_id:
                        continue
                    tmp_sheet = self.sheet_objs[sid]
                    if check_common_elements(cur_face.Eids, tmp_sheet.parallel_eids):
                        continue
                    neighbor_sid.add(sid)
                    sheet.neighbor_sheet_ids.add(sid)
                    tmp_sheet.neighbor_sheet_ids.add(sheet.sheet_id)
            wall_neighbor_sheet_ids[i] = neighbor_sid
        sheet.separate_wall_neighbor_sids = wall_neighbor_sheet_ids

    def build_sheet_neighbor_connectivity(self):
        print("build sheet neighbor connectivity ...")
        for sheet in tqdm(self.sheet_objs):
            self.find_sheet_neighbor_sheet(sheet)
        self.assign_sheet_colors()
        print("done ...")

    def get_smaller_one_side_sheets(self, sheet_id):
        separate_neighbor_sids = self.sheet_objs[sheet_id].separate_wall_neighbor_sids
        neighbor_cell_number = []
        valid_sid = []
        for wall_id in separate_neighbor_sids:
            nsids = separate_neighbor_sids[wall_id]
            if len(nsids) == 0:
                continue
            valid_sid.append(wall_id)
            cell_number = 0
            for sid in nsids:
                cell_number += len(self.sheet_objs[sid].hexa_cids)
            neighbor_cell_number.append(cell_number)
        if not neighbor_cell_number:
            return -1
        idx = np.argmin(neighbor_cell_number)
        idx = valid_sid[idx]
        return separate_neighbor_sids[idx]

    # -------------- sheet relationship building process done -------------------
    # ToDo:
    # each parallel edge should have its own parallel edge link of lower level
    def get_sheet_background_parallel_eids(self, sheet_id):
        sheet = self.sheet_objs[sheet_id]
        current_mesh = self.mesh
        cur_parallel_eids = []
        for eid in sheet.parallel_eids:
            if current_mesh.mesh_level > 0:
                cur_parallel_eids.extend(current_mesh.Edges[eid].es_link)
            else:
                cur_parallel_eids.append(eid)
        while current_mesh.mesh_level > 0:
            current_mesh = current_mesh.input_mesh
            next_level_flag = current_mesh.mesh_level
            next_parallel_eids = set()
            visited_eids = set()
            while cur_parallel_eids:
                eid = cur_parallel_eids.pop()
                if eid in visited_eids:
                    continue
                visited_eids.add(eid)
                edge = current_mesh.Edges[eid]
                cur_parallel_eids.extend(edge.parallels)
                if next_level_flag:
                    next_parallel_eids.update(edge.es_link)
            if next_level_flag:
                cur_parallel_eids = list(next_parallel_eids)
            else:
                cur_parallel_eids = list(visited_eids)
        return current_mesh, cur_parallel_eids

    # invalid matching is possible
    def create_matching_map(self):
        for s in self.sheet_objs:
            eids = s.parallel_eids
            for e in eids:
                vids = self.mesh.Edges[e].Vids
                for vid in vids:
                    if vid not in s.vert_matching_map:
                        s.vert_matching_map[vid] = set()
                    s.vert_matching_map[vid].update(vids)
                    # remove self from set
                    s.vert_matching_map[vid].remove(vid)
            self.verify_sheet_is_prefect(s.sheet_id)
            self.find_mul_matching_cells(s.sheet_id)
            self.find_self_intersecting_cells(s.sheet_id)
            self.find_perfect_edge(s.sheet_id)

    def find_perfect_edge(self, sheet_id):
        sheet = self.sheet_objs[sheet_id]
        for eid in sheet.parallel_eids:
            vids = self.mesh.Edges[eid].Vids
            is_perfect = True
            for vid in vids:
                if len(sheet.vert_matching_map[vid]) != 1:
                    is_perfect = False
                    break
            if is_perfect:
                sheet.prefect_parallel_edge.add(eid)

    def verify_sheet_is_prefect(self, sheet_id):
        sheet = self.sheet_objs[sheet_id]
        perfect_flag = True
        for k, v in sheet.vert_matching_map.items():
            if len(v) != 1:
                perfect_flag = False
                break
        sheet.is_perfect = perfect_flag

    def find_mul_matching_cells(self, sheet_id):
        sheet = self.sheet_objs[sheet_id]
        if sheet.is_perfect:
            return None
        mul_matching_cids = set()
        # only hex cells
        # print(sheet.hexa_cids)
        for cid in sheet.cells:
            cell = self.mesh.Cells[cid]
            if not cell.is_hexa():
                continue
            for v in cell.Vids:
                if len(sheet.vert_matching_map[v]) == 1:
                    continue
                mul_matching_cids.add(cid)
                break
        sheet.mul_matching_cids = mul_matching_cids
        return sheet.mul_matching_cids

    def find_self_intersecting_cells(self, sheet_id):
        sheet = self.sheet_objs[sheet_id]
        if sheet.is_perfect:
            return None
        self_intersecting_cids = set()
        for cid in sheet.cells:
            cell = self.mesh.Cells[cid]
            if not cell.is_hexa():
                continue
            in_sheet_eids = check_common_elements(sheet.parallel_eids, cell.Eids)
            for eid in in_sheet_eids:
                edge = self.mesh.Edges[eid]
                neids = check_common_elements(edge.neighboring_Eids, in_sheet_eids)
                if len(neids) > 0:
                    self_intersecting_cids.add(cid)
                    break
                # if cid in self_intersecting_cids:
                #     break
        sheet.self_intersecting_cids = self_intersecting_cids
        # print(list(self_intersecting_cids))

    def set_max_sheet_neighbor_num(self):
        max_sheet_neighbor_number = 0
        for sheet in self.sheet_objs:
            for wall_id in sheet.separate_wall_neighbor_sids:
                wall_sids = sheet.separate_wall_neighbor_sids[wall_id]
                v = len(wall_sids)
                if v > max_sheet_neighbor_number:
                    max_sheet_neighbor_number = v
        self.max_neighbor_num = max_sheet_neighbor_number
        self.color_num = 2 + max_sheet_neighbor_number

    def assign_sheet_colors(self):
        self.set_max_sheet_neighbor_num()
        for sheet in self.sheet_objs:
            neighbor_colors = set()
            for sid in sheet.neighbor_sheet_ids:
                if self.sheet_objs[sid].color_id:
                    neighbor_colors.add(self.sheet_objs[sid].color_id)
            candidate_colors = set()
            for i in range(self.color_num):
                # keep rad color reserved
                if i == 0:
                    continue
                if i in neighbor_colors:
                    continue
                candidate_colors.add(i)
            if not candidate_colors:
                self.color_num += 1
                candidate_colors.add(self.color_num)
            # print(list(candidate_colors))
            color_list = list(candidate_colors)
            picked_color = random.choice(color_list)
            # picked_color = list(candidate_colors)[int(len(candidate_colors) / 2)]
            sheet.color_id = picked_color

    def element_check(self):
        # check section
        sheet = self.sheet_objs[4]
        cell = self.mesh.Cells[1026]
        print(f"cell eids : {cell.Eids}")
        print(f"vert edges : {self.mesh.Vertices[826].neighboring_Eids}")
        for eid in self.mesh.Vertices[826].neighboring_Eids:
            edge = self.mesh.Edges[eid]
            print(f"edge {eid} parallel edges")
            print(edge.parallels)
            for ep in edge.parallels:
                eee = self.mesh.Edges[ep]
                print(f"edge {eid} parallel edge {ep}'s parallel edges")
                print(eee.parallels)
            print("---")
        print(f"cell in sheet eids : {check_common_elements(sheet.parallel_eids, cell.Eids)}")

        # check section done

    def get_largest_sheet_ids(self):
        max_sids = set()
        max_cell_size = 0
        for i, se in enumerate(self.sheet_objs):
            sheet_size = len(se.cells)
            if sheet_size > max_cell_size:
                max_sids = set()
                max_cell_size = sheet_size
            if sheet_size == max_cell_size:
                max_sids.add(i)
        return max_sids

    def get_size_ranked_imperfect_sheet_ids(self):
        imperfect_sids = self.get_imperfect_sheet_ids()
        sheet_size_dict = {i: 0 for i in imperfect_sids}
        for i in imperfect_sids:
            se = self.sheet_objs[i]
            sheet_size = len(se.cells)
            sheet_size_dict[i] = sheet_size
        return sorted(sheet_size_dict)

    def get_imperfect_sheet_ids(self):
        imperfect_sids = [s.sheet_id for s in self.sheet_objs if not s.is_perfect]
        return imperfect_sids

    def get_perfect_sheet_ids(self):
        perfect_sids = [s.sheet_id for s in self.sheet_objs if s.is_perfect]
        return perfect_sids

    def get_size_ranked_sheet_ids(self):
        sheet_size_dict = {i: 0 for i in range(len(self.sheet_objs))}
        for i in sheet_size_dict:
            se = self.sheet_objs[i]
            sheet_size = len(se.cells)
            sheet_size_dict[i] = sheet_size
        return sorted(sheet_size_dict)

    def mesh_cover_max_sheet_set(self):
        grouped_sheet_ids = set()
        visited_cells = set()
        hex_cells = set(self.mesh.get_hexa_cell_id_list())
        unvisited_hex_cell = hex_cells
        ranked_sids = self.get_size_ranked_sheet_ids()
        # remove self-intersecting sheets
        self_intersecting_sids = set()
        for sid in ranked_sids:
            sheet = self.sheet_objs[sid]
            if sheet.self_intersecting_cids:
                self_intersecting_sids.add(sid)
        # remove
        for sid in self_intersecting_sids:
            ranked_sids.remove(sid)
        candidate_sids = [ranked_sids[0]]
        while candidate_sids:
            sid = candidate_sids.pop()
            sheet = self.sheet_objs[sid]
            ranked_sids.remove(sid)
            grouped_sheet_ids.add(sid)
            visited_cells.update(sheet.cells)
            unvisited_hex_cell = unvisited_hex_cell.difference(visited_cells)
            if not unvisited_hex_cell:
                break
            if not ranked_sids:
                break
            # find next sheet
            overlap_cell_number = []
            for nsid in ranked_sids:
                n_sheet = self.sheet_objs[nsid]
                overlap_cell_n = len(check_common_elements(n_sheet.cells, visited_cells))
                overlap_cell_number.append(overlap_cell_n)
            # next number
            if 0 in overlap_cell_number:
                indexes = overlap_cell_number.index(0)
                candidate_sids.append(ranked_sids[indexes])
            else:
                min_pos = overlap_cell_number.index(min(overlap_cell_number))
                candidate_sids.append(ranked_sids[min_pos])
        return grouped_sheet_ids











