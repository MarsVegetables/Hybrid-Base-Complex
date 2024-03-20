class ParallelSheetAdjacencyInfo:
    def __init__(self, sid):
        self.parallel_sid = sid
        self.parallel_candidate_non_hexa_cids = {}
        self.neighbor_sids = {}  # {sid : corresponding non hexa cids}
        self.cell_merge_groups = []

    def add_neighbor_sid_and_correspond_non_hexa_cids(self, neighbor_sid, non_hexa_cids):
        if neighbor_sid in self.neighbor_sids:
            self.neighbor_sids[neighbor_sid].update(non_hexa_cids)
        else:
            self.neighbor_sids[neighbor_sid] = set(non_hexa_cids)


# sheet contain non hex or hex cells
class HybridSheet:
    def __init__(self, sheet_id):
        self.sheet_id = sheet_id
        self.parallel_eids = []
        # cells include non hex and hex cell ids
        self.cells = []
        self.boundary_fids = set()  # sheet boundary fids
        self.wall_fids = []  # the faces that does not contain parallel edge
        self.separate_wall_fids = []  # each column is a set of fid in same wall
        self.separate_wall_neighbor_sids = []
        # inner face : faces that inside sheet
        self.sheet_inner_fids = set()
        # inner fids patch may include boundary face if they have connection
        self.separate_inner_fids_patch = []
        self.hexa_cids = set()  # only contain hexa
        self.neighbor_sheet_ids = set()
        self.full_neighbor_sheet_ids = set()  # all member are neighboring with a sheet in this list
        self.partial_neighbor_sheet_ids = set()  # part of members are neighboring with a sheet in this list
        self.non_hexa_cids = set()
        # vid : connect vids
        self.vert_matching_map = {}
        # perfect parallel edge
        self.prefect_parallel_edge = set()
        self.mul_matching_cids = set()
        self.self_intersecting_cids = set()
        self.is_perfect = False
        self.color_id = None

    def color(self):
        return self.color_id

    def get_unmatched_vids(self):
        vids = set()
        for vid in self.vert_matching_map:
            vm = self.vert_matching_map[vid]
            if len(vm) != 1:
                vids.add(vid)
        return vids
