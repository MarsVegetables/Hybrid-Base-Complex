'''
Date: 2021-08-08 22:06:47
LastEditors: Lei Si
LastEditTime: 2023-02-20 14:44:02
'''

Interior_Regular_E = 4
Boundary_Regular_E = 2

ELEMENT_TYPE_LIST = ["Base", # 0
                     "Vertex", # 1
                     "Edge", # 2
                     "Face", # 3
                     "Cell", # 4
                     "Tet"  # 5
                     "Hexa"  # 6
                     ]


class BlenderProperties:
    def __init__(self):
        # FIFD
        # for a component, multiple background elements need to be displayed one by one
        # background element start time = component key frame time + index in element_display_order
        self.element_display_order = []
        # key_frame is the time of display
        self.key_frame = 0
        self.color = 0  # keep use newest color
        self.old_color = -1  # update when merge to neighbor sheets


class SuperProperties:
    def __init__(self):
        # ------------------------------------------
        # super element attributes
        # ------------------------------------------
        # the basic member of this element is saved in basic_elements
        self.basic_elements = []
        # for a super element
        # connect elements type
        # for edge is vert
        # for face is edge
        # for cell is face
        self.connect_elements = []
        # if a element is basic element
        # use super id to find corresponding super element
        # a basic element can has multiple super elements
        self.super_ids = set()

        # true if the element is super component element
        self.is_super = False

    def add_basic_element_ids(self, element_ids):
        for i in element_ids:
            if i in self.basic_elements:
                continue
            self.basic_elements.append(i)

    # connect element is n-1 cell
    # for example,
    # a connect element of a face (2 cell) must be a vertex (1 cell)
    def add_connect_element_ids(self, element_ids):
        for i in element_ids:
            if i in self.connect_elements:
                continue
            self.connect_elements.append(i)


class ComponentProperties:
    def __init__(self):
        # base complex component attr
        # self.component_id = -1
        # for t-mesh, mesh element may include in more than one component elements
        self.component_ids = set()
        # number of element
        self.resolution = 1
        # True is the component need to be refined.
        self.need_refine = False

    def is_component(self):
        return len(self.component_ids) != 0


class ElementBase(BlenderProperties, SuperProperties, ComponentProperties):
    def __init__(self, element_id, element_type=0):
        BlenderProperties.__init__(self)
        SuperProperties.__init__(self)
        ComponentProperties.__init__(self)
        self.element_type = element_type
        self.element_id = element_id
        self.visited = False
        self.removed = 0

        # full neighbor
        self.neighboring_Vids = []
        self.neighboring_Eids = []
        self.neighboring_Fids = []
        self.neighboring_Cids = []

        # all neighbor ids
        self.neighboring_Vids_a = []
        self.neighboring_Eids_a = []
        self.neighboring_Fids_a = []
        self.neighboring_Cids_a = []

        # partial neighbor
        self.neighboring_Vids_p = []
        self.neighboring_Eids_p = []
        self.neighboring_Fids_p = []
        self.neighboring_Cids_p = []

        self.is_boundary = False
        self.is_sharp = False
        self.is_singular = False
        self.is_non_manifold = False

        # True if
        # cell is hex
        # or
        # face, edge, vertex
        # are only neighboring with hex cells
        self.hex_flag = False

    def print_id_type(self):
        print("type : " + ELEMENT_TYPE_LIST[self.element_type])
        print("id : " + str(self.element_id))

    def print_self_data(self):
        self.print_id_type()

    # combine full neighbor and partial neighbor list
    def get_all_neighbor(self):
        pass

    def print_neighboring(self, id_flag=True):
        if id_flag:
            self.print_id_type()

        print("neighboring_Vids --> " + str(self.neighboring_Vids))
        print("neighboring_Eids --> " + str(self.neighboring_Eids))
        print("neighboring_Fids --> " + str(self.neighboring_Fids))
        print("neighboring_Cids --> " + str(self.neighboring_Cids))

        self.print_partial_neighbor()

    def print_partial_neighbor(self):
        print("partial neighboring Vids --> " + str(self.neighboring_Vids_p))
        print("partial neighboring Eids --> " + str(self.neighboring_Eids_p))
        print("partial neighboring Fids --> " + str(self.neighboring_Fids_p))
        print("partial neighboring Cids --> " + str(self.neighboring_Cids_p))

    def print_all(self):
        self.print_self_data()
        self.print_neighboring(False)

    def init_neighbors(self):
        # init or reset
        self.neighboring_Vids = []
        self.neighboring_Eids = []
        self.neighboring_Fids = []
        self.neighboring_Cids = []

        # all neighbor ids
        self.neighboring_Vids_a = []
        self.neighboring_Eids_a = []
        self.neighboring_Fids_a = []
        self.neighboring_Cids_a = []

        # partial neighbor
        self.neighboring_Vids_p = []
        self.neighboring_Eids_p = []
        self.neighboring_Fids_p = []
        self.neighboring_Cids_p = []


#   
#
class Vertex(ElementBase):
    def __init__(self, element_id, coordinate, element_type = 1):
        super(Vertex, self).__init__(element_id, element_type)
        self.x = coordinate[0]
        self.y = coordinate[1]
        self.z = coordinate[2]

    def init_svid_fvid(self):
        # base complex attribute
        self.svid = -1  # singular vert id
        self.fvid = -1  # frame vert id

    def print_self_data(self):
        super().print_self_data()
        print("coordinate : " + str([self.x, self.y, self.z]))

    def xyz(self):
        return [self.x, self.y, self.z]
    
    def add_neighboring_vids(self, edges):
        for edge in edges:
            for nvid in edge:
                if nvid not in self.neighboring_Vids and nvid != self.element_id:
                    self.neighboring_Vids.append(nvid)

    def add_neighboring_eids(self, eids):
        for eid in eids:
            if eid not in self.neighboring_Eids:
                self.neighboring_Eids.append(eid)
    
    def add_neighboring_fids(self, fids):
        for fid in fids:
            if fid not in self.neighboring_Fids:
                self.neighboring_Fids.append(fid)

    def add_neighboring_cids(self, cids):
        for cid in cids:
            if cid not in self.neighboring_Cids:
                self.neighboring_Cids.append(cid)

    def update_xyz(self, new_coordinate):
        self.x = new_coordinate[0]
        self.y = new_coordinate[1]
        self.z = new_coordinate[2]


# edge type
#
class Edge(ElementBase):
    def __init__(self, element_id, vids, element_type = 2):
        super(Edge, self).__init__(element_id, element_type)
        self.Vids = vids
        self.parallels = []
        # base complex attribute
        self.hex_edge = False
        self.seid = -1 # singular edge id
        self.T_edge = False

    def print_self_data(self):
        super().print_self_data()
        print("Vids : " + str(self.Vids))

    def add_neighboring_vids(self, vids):
        for vid in vids:
            if vid not in self.Vids:
                if vid not in self.neighboring_Vids:
                    self.neighboring_Vids.append(vid)

    def add_neighboring_eids(self, eids):
        for eid in eids:
            if eid != self.element_id and eid not in self.neighboring_Eids:
                self.neighboring_Eids.append(eid)

    def add_neighboring_fids(self, fids):
        for fid in fids:
            if fid not in self.neighboring_Fids:
                self.neighboring_Fids.append(fid)
    
    def add_neighboring_cids(self, cids):
        for cid in cids:
            if cid not in self.neighboring_Cids:
                self.neighboring_Cids.append(cid)
    
    def add_parallels(self, elementId):
        # this function does not varify the correctness.
        if elementId == self.element_id:
            return
        if elementId not in self.parallels:
            self.parallels.append(elementId)
    
    def is_parallel(self, eid):
        if eid in self.parallels:
            return True
        return False


#
#
#
class Face(Edge):
    def __init__(self, element_id, vids, eids, element_type = 3):
        super(Face, self).__init__( element_id, vids, element_type)
        self.Eids = eids
        self.oppositeEdgeFacePairs = {eid: [] for eid in self.Eids}
        # self.all_hexa_neighbor = 0
        # true if the face connect two cell to build a super cell
        self.isConnectFace = False
        # sheet information
        self.neighboring_sheet_ids = set()

    def update_eids(self, eids):
        self.Eids = eids
        self.oppositeEdgeFacePairs = {eid: [] for eid in self.Eids}

    def print_self_data(self):
        super().print_self_data()
        print("Eids : " + str(self.Eids))
        print("Opposite Faces : " + str(self.oppositeEdgeFacePairs))
        print("parallel Faces : " + str(self.parallels))

    def add_neighboring_eids(self, eids):
        for eid in eids:
            if eid not in self.Eids:
                if eid not in self.neighboring_Eids:
                    self.neighboring_Eids.append(eid)

    def add_neighboring_fids(self, fids):
        for fid in fids:
            if fid != self.element_id and fid not in self.neighboring_Fids:
                self.neighboring_Fids.append(fid)     
    
    def add_opposite_edge_face_pairs(self, eid, fid):
        # this function does not varify the correctness.
        if eid not in self.Eids:
            print("Eid : " + str(eid) + " is not included in face : " + str(self.element_id))
            return
        if fid == self.element_id:
            return
        if fid not in self.oppositeEdgeFacePairs[eid]:
            self.oppositeEdgeFacePairs[eid].append(fid)
    
    def add_parallel_faces(self, fid):
        # this function does not varify the correctness.
        self.add_parallels(fid)

    def is_quad(self):
        quad_flag = True
        if len(self.Eids) != 4:
            quad_flag = False
        if len(self.Vids) != 4:
            quad_flag = False
        return quad_flag


#
#
class Cell(Face):
    def __init__(self, element_id, vids, eids, fids, fOrientations = [], cellOrientation = 0, element_type = 4):
        super(Cell, self).__init__(element_id, vids, eids, element_type)
        self.Fids = fids
        self.faceOrientations = fOrientations
        self.cellOrientation = cellOrientation

    def print_self_data(self):
        super().print_self_data()
        print("Fids : " + str(self.Fids))

    def add_neighboring_fids(self, fids):
        for fid in fids:
            if fid not in self.Fids:
                if fid not in self.neighboring_Fids:
                    self.neighboring_Fids.append(fid)

    def add_neighboring_cids(self, cids):
        for cid in cids:
            if cid != self.element_id and cid not in self.neighboring_Cids:
                self.neighboring_Cids.append(cid)
                
    def corner_vids(self, vid):
        pass
    
    def is_hexa(self):
        if self.element_type == 6:
            return True
        return False
    
    def get_hybrid_cell_format(self):
        if len(self.faceOrientations) == 0:
            self.faceOrientations = [0] * len(self.Fids)
        return [self.Fids, self.faceOrientations, self.cellOrientation]