from mesh_structure.mesh_elements import ElementBase, Cell, Face, Edge, Vertex


class SingularElementProperties:
    def __init__(self):
        self.is_fake = False


class SingularityVert(ElementBase, SingularElementProperties):
    def __init__(self, element_id, mesh_vid, element_type=1):
        ElementBase.__init__(self, element_id, element_type)
        SingularElementProperties.__init__(self)
        self.mesh_vid = mesh_vid

    def init_attrs(self):
        self.is_fake = False


class SingularityEdge(ElementBase, SingularElementProperties):
    def __init__(self, element_id, element_type=2):
        ElementBase.__init__(self, element_id, element_type)
        SingularElementProperties.__init__(self)
        self.init_attrs()

    def init_attrs(self):
        self.startEndVids = [] # singularity vid
        self.es_link = []  # hex e id
        self.vs_link = []  # hex v id
        self.is_boundary = False
        self.is_fake = False
        self.is_circle = False
        self.componentEids_link = [] # base-complex es, one or more base-complex edges compose one singular edge
        self.separatedFacePatchIds = set()
        # self.neighborComponentCidsGroups = []


class ComponentCell(Cell):
    def __init__(self, element_id, cids = []):
        Cell.__init__(self, element_id=element_id, vids=[], eids=[], fids=[])
        self.cids_patch = cids
        self.face_eids = {}


class ComponentFace(Face):
    def __init__(self, element_id, fids=[], is_boundary=False):
        Face.__init__(self, element_id, vids=[], eids=[])
        # mesh fids
        self.fids_patch = fids
        self.is_boundary = is_boundary
        self.neighboring_sheet_ids = set()  # store base complex sheet ids
        self.ranked_predicted_valence_eids = {}


class ComponentEdge(Edge, SingularityEdge):
    def __init__(self, element_id, es_link=[], is_boundary=False):
        Edge.__init__(self, element_id, vids=[])
        SingularityEdge.__init__(self, element_id)
        self.es_link = es_link
        self.vs_link = []
        self.is_boundary = is_boundary


class ComponentVert(Vertex, SingularityVert):
    def __init__(self, element_id, mesh_vid=-1, xyz=[-1,-1,-1], is_boundary=False):
        Vertex.__init__(self, element_id, coordinate=xyz)
        SingularityVert.__init__(self, element_id, mesh_vid)
        self.is_boundary = is_boundary


class candidateNonHexaCell:
    def __init__(self, element_id, cid):
        self.element_id = element_id
        self.cid = cid
        self.merge_eids = []
        self.candidate_fids = []
        self.candidate_neighboring_hexa_cids = [] # neighbor with the candidate fids
        self.parallel_sids = {} # sid : parallel eids
        self.neighboring_sids = {} # sid : hexa cids

    def __str__(self):
        info = f"""
        {'-'*40}
        # id : {self.element_id}
        # cid : {self.cid}
        # merge eids : {self.merge_eids}
        # candidate fids : {self.candidate_fids}
        # candidate neighboring hexa cids : {self.candidate_neighboring_hexa_cids}
        # parallel sid : {self.parallel_sids}
        # neighbor sids : {self.neighboring_sids}
        {'-'*40}
        """
        return info


class parallelSheetAdjacencyInfo:
    def __init__(self, sid):
        self.parallel_sid = sid
        self.parallel_candidate_non_hexa_cids = {}
        self.neighbor_sids = {} # {sid : corresponding non hexa cids}
        self.cell_merge_groups = []

    def addNeighborSidAndCorrespondNonHexaCids(self, neighbor_sid, non_hexa_cids):
        if neighbor_sid in self.neighbor_sids:
            self.neighbor_sids[neighbor_sid].update(non_hexa_cids)
        else:
            self.neighbor_sids[neighbor_sid] = set(non_hexa_cids)


# base complex sheet only contain hexa component cell
# a hexa component cell may include few non hexa elements
class hybridBaseComplexSheet:
    def __init__(self, sheet_id):
        self.sheet_id = sheet_id
        self.parallelEids = []
        self.boundary_component_fids = set() # sheet is_boundary component fids
        self.wall_component_fids = [] # the faces that does not contain parallel edge
        self.separate_wall_component_fids = [] # each column is a set of fid in same wall
        self.separate_wall_neighbor_sids = []
        self.hexa_component_cids = [] # onlt contain hexa base complex component
        self.neighbor_sheet_ids = set()
        self.full_neighbor_sheet_ids = set() # all member are neighboring with a sheet in this list
        self.partial_neighbor_sheet_ids = set() # part of members are neighboring with a sheet in this list
        self.neighbor_non_hexa_component_cids = []



