from mesh_structure.mesh_common import check_common_elements


class RegionCellTracer:
    def __init__(self, mesh, parallel_eids, region_cids):
        self.mesh = mesh
        self.parallel_eids = parallel_eids
        self.region_cids = region_cids
        self.extend_cids = {}

    # extend all neighboring non-hex cells
    # may go to other layer
    def extend_non_hex_cells(self):
        self.extend_cids = set()
        eid_stack = list(self.parallel_eids)
        visited_eids = set()
        while eid_stack:
            eid = eid_stack.pop()
            if eid in visited_eids:
                continue
            visited_eids.add(eid)
            edge = self.mesh.Edges[eid]
            for cid in edge.neighboring_Cids:
                if cid in self.region_cids:
                    continue
                cell = self.mesh.Cells[cid]
                if cell.is_hexa():
                    continue
                self.extend_cids.add(cid)
                eid_stack.extend(cell.Eids)
        return self.extend_cids

    def extend_non_hex_cells_one_ring(self):
        self.extend_cids = set()
        eid_stack = list(self.parallel_eids)
        visited_eids = set()
        while eid_stack:
            eid = eid_stack.pop()
            if eid in visited_eids:
                continue
            visited_eids.add(eid)
            edge = self.mesh.Edges[eid]
            for cid in edge.neighboring_Cids:
                if cid in self.region_cids:
                    continue
                cell = self.mesh.Cells[cid]
                if cell.is_hexa():
                    continue
                self.extend_cids.add(cid)
        return self.extend_cids


