from structure_analysis.sheet_extraction import SheetExtraction
from mesh_structure.mesh_common import check_common_elements, reduce_id_groups


class SheetIdentification(SheetExtraction):
    def __init__(self, mesh):
        SheetExtraction.__init__(self, mesh)
        self.vert_weight = {i: 0 for i in self.mesh.Vertices}

    def get_sheet_imperfect_element_ids(self, sid):
        sheet = self.sheet_objs[sid]
        if sheet.is_perfect:
            return set(), set(), set()

        self_intersecting_cids = set()
        self_parallel_eids = set()
        unmatched_vids = set()

        for vid, value in sheet.vert_matching_map.items():
            if len(value) != 1:
                self.vert_weight[vid] += 1
                unmatched_vids.add(vid)

        # detect self intersecting
        for cid in sheet.cells:
            in_sheet_eids = check_common_elements(self.mesh.Cells[cid].Eids, sheet.parallel_eids)
            if len(in_sheet_eids) != 4:
                self_intersecting_cids.add(cid)

        # detect self parallel
        # also complete the self intersecting (redundent)
        for eid in sheet.parallel_eids:
            edge = self.mesh.Edges[eid]
            tmp = check_common_elements(edge.neighboring_Eids, sheet.parallel_eids)
            for neid in tmp:
                n_edge = self.mesh.Edges[neid]
                common_cids = check_common_elements(edge.neighboring_Cids, n_edge.neighboring_Cids)
                if common_cids:
                    self_intersecting_cids.update(check_common_elements(common_cids, sheet.cells))
                else:
                    self_parallel_eids.add(neid)
                    self_parallel_eids.add(eid)
        return self_intersecting_cids, self_parallel_eids, unmatched_vids

    def get_two_sheet_similarity(self, sheet_id_1, sheet_id_2, neighbor_ratio):
        sheet_1 = self.sheet_objs[sheet_id_1]
        sheet_2 = self.sheet_objs[sheet_id_2]
        # shared fids
        shared_fids = check_common_elements(sheet_1.wall_fids, sheet_2.wall_fids)
        # ratio

        fid_count = 0
        for wf in sheet_1.separate_wall_fids:
            if not check_common_elements(wf, shared_fids):
                continue
            fid_count += len(wf)
        s1_max = len(shared_fids) / fid_count

        fid_count = 0
        for wf in sheet_2.separate_wall_fids:
            if not check_common_elements(wf, shared_fids):
                continue
            fid_count += len(wf)
        s2_max = len(shared_fids) / fid_count
        r = min(s1_max, s2_max)
        if r >= neighbor_ratio:
            return r
        return 0

    def sheet_grouping(self, neighbor_ratio):
        sheet_groups = []
        for i, sheet_obj in enumerate(self.sheet_objs):
            g = {i}
            nsids = sheet_obj.neighbor_sheet_ids
            for nsid in nsids:
                r = self.get_two_sheet_similarity(i, nsid, neighbor_ratio)
                if r > 0:
                    g.add(nsid)
            sheet_groups.append(list(g))
        sheet_groups = reduce_id_groups(sheet_groups)
        return sheet_groups




