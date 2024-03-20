import random

import numpy as np

from .mesh_elements import Vertex, Edge, Face, Cell, Interior_Regular_E, Boundary_Regular_E
from .mesh_common import (face_list_to_edge_list, check_common_elements,
                          reduce_id_groups, remove_duplicate_faces_and_remapping,
                          find_duplicate_elements)
from tqdm import tqdm
import pathlib


class BasicMesh:
    # ToDo: split this large class into multiple shorter classes
    def __init__(self,  file_path, mesh_type) -> None:
        # file path initialize
        self.file_path = ""
        self.file_format = ""
        self.meshType = ""
        self.structure_type = "basic"
        # mesh structure initialize
        self.Vertices = {}
        self.Edges = {}
        self.Faces = {}
        self.Cells = {}
        self.max_poly_face_num = 0
        # color number determined by max number of neighbor cells of a cell
        self.color_num = 0
        # edge color
        self.edge_color_num = 0

        self.illnessFaces = set()
        self.illnessCells = set()
        self.illnessEdges = set()
        self.multi_opposite_fids = set()  # if a face contain more than one opposite face
        self.include_components = False
        self.hexa_singular_eids = set()
        self.non_hexa_eids = set()

        # assign path and mesh type
        self.update_file_path(file_path)
        self.update_mesh_type(mesh_type)

        # ToDo : create a matrix (graph) during load data may can improve performance

    # buildMeshFromVFC
    def build_mesh_from_vfc(self, points, faces, cells):
        # require:
        #  vertice xyz list
        #  face vids list
        #  cell fid fos co list
        # must do not sort input face vids.
        cleaned_faces, cleaned_cells = remove_duplicate_faces_and_remapping(faces, cells)
        self.init_attrs()
        print("init mesh elements ...")
        self.init_vertices(points)
        edges, faces_eids = face_list_to_edge_list(cleaned_faces)
        self.init_edges(edges)
        self.init_faces(cleaned_faces, faces_eids)
        self.init_cells_with_fid_and_orientation(cleaned_cells)
        print("done ...")
        self.build()
        # self.show_statistic()

    def init_attrs(self):
        self.Vertices = {}
        self.Edges = {}
        self.Faces = {}
        self.Cells = {}
        self.max_poly_face_num = 0
        self.color_num = 0
        self.edge_color_num = 0
        self.illnessFaces = set()
        self.illnessCells = set()
        self.illnessEdges = set()
        self.multi_opposite_fids = set()  # if a face contain more than one opposite face
        self.include_components = False
        self.hexa_singular_eids = set()
        self.non_hexa_eids = set()

    def update_file_path(self, f_path):
        self.file_path = f_path
        file_ext = pathlib.Path(f_path).suffix
        self.file_format = file_ext

    def update_mesh_type(self, meshType):
        self.meshType = meshType

    # -------------- init element -------------------
    def update_vertices(self, vid, vclass):
        self.Vertices.update({vid : vclass})

    def add_vertex_element(self, vid, coord):
        if vid not in self.Vertices.keys():
            vc = Vertex(vid, coord)
            self.update_vertices(vid, vc)

    def update_vertex_element(self, vid, coord):
        vc = Vertex(vid, coord)
        self.update_vertices(vid, vc)

    def init_vertices(self, points):
        # input : vertex coordination list
        #
        # add vertex element
        # vertex element connectivies will be complete when this function finish.
        #
        print("init Vertex ...")
        for vid, coord in enumerate(points):
            self.add_vertex_element(vid, coord)
        print(f"init {len(points)} Vertices ... done")

    def update_edges(self, eid, eclass):
        self.Edges.update({eid: eclass})

    def add_edge_element(self, eid, vids):
        if eid not in self.Edges.keys():
            ec = Edge(eid, list(vids))
            self.update_edges(eid, ec)

    def update_edge_element(self, eid, vids):
        ec = Edge(eid, list(vids))
        self.update_edges(eid, ec)

    def init_edges(self, edges):
        print("init Edges ...")
        for i, e_vids in enumerate(edges):
            self.add_edge_element(i, e_vids)
        print(f"init {len(edges)} Edges ... done")

    def is_face_duplicate(self, fid, face_vids):
        for existing_fid in self.Faces:
            face = self.Faces[existing_fid]
            vids = set(face.Vids)
            if set(face_vids) == set(vids):
                return existing_fid
        return fid

    def update_faces(self, fid, fclass):
        self.Faces.update({fid : fclass})

    def add_face_element(self, fid, vids, eids):
        if fid not in self.Faces.keys():
            fc = Face(fid, list(vids), eids)
            self.update_faces(fid, fc)

    def update_face_element(self, fid, vids, eids):
        fc = Face(fid, list(vids), eids)
        self.update_faces(fid, fc)

    def init_faces(self, faces, faces_eids):
        print("init Faces ...")
        nf = len(faces)
        if nf != len(faces_eids):
            print("face number is not match face faces_eids number : " + str(nf) + " -- " + str(faces_eids))
        for fid in range(nf):
            self.add_face_element(fid, faces[fid], faces_eids[fid])
        print(f"init {nf} Faces ... done")

    def update_cells(self, cid, cclass):
        self.Cells.update({cid : cclass})

    def add_cell_element(self, cid, vids, eids, fids, fos = [], co = 0):
        # fos is face orientation in this cell
        # co is cell orientation
        if cid not in self.Cells.keys():
            cc = Cell(cid, vids, eids, fids, fos, co)
            self.update_cells(cid, cc)

    def update_cell_element(self, cid, vids, eids, fids, fos = [], co = 0):
        # fos is face orientation in this cell
        # co is cell orientation
        cc = Cell(cid, vids, eids, fids, fos, co)
        self.update_cells(cid, cc)

    def init_cells_with_fid_and_orientation(self, cells):
        # cell list structure
        #   [[face id list], 
        #    [face orientation list], 
        #    cell orientation]
        # if value of orientation is 1, the face or cell need to change the normal direction.
        print("init Cells ...")
        for cid, cell in enumerate(cells):
            fids = cell[0]
            fos = cell[1]
            co = cell[2]
            self.add_cell_element(cid, [], [], fids, fos, co)
            if len(fids) > self.max_poly_face_num:
                self.max_poly_face_num = len(fids)
        print(f"init {len(cells)} Cells ... done")

    def init_vert_svid_fvid(self):
        for vid, vert in self.Vertices.items():
            vert.init_svid_fvid()

    def clean_visited(self):
        for vid in self.Vertices.keys():
            self.Vertices[vid].visited = 0
        for fid in self.Faces.keys():
            self.Faces[fid].visited = 0
        for eid in self.Edges.keys():
            self.Edges[eid].visited = 0
        for cid in self.Cells.keys():
            self.Cells[cid].visited = 0

    def clean_removed(self):
        for vid in self.Vertices.keys():
            self.Vertices[vid].removed = 0
        for fid in self.Faces.keys():
            self.Faces[fid].removed = 0
        for eid in self.Edges.keys():
            self.Edges[eid].removed = 0
        for cid in self.Cells.keys():
            self.Cells[cid].removed = 0

    def clean_neighbors(self):
        print("cleanning element neighbors ...")
        for vid in self.Vertices.keys():
            self.Vertices[vid].init_neighbors()
        for fid in self.Faces.keys():
            self.Faces[fid].init_neighbors()
        for eid in self.Edges.keys():
            self.Edges[eid].init_neighbors()
        for cid in self.Cells.keys():
            self.Cells[cid].init_neighbors()
        print("cleanning element neighbors ... done")

    # -------------- init element done -------------------

    # -------------- build connectivity -------------------
    def build_vertex_connectivity(self):
        # start from cell
        # verify hexa cell at the same time.
        print("building vertex connectivity ...")
        for cid in tqdm(self.Cells):
            if self.Cells[cid].removed:
                continue
            self.build_vertex_connectivity_cid(cid)
        print("done ...")

    def build_vertex_connectivity_cid(self, cid):
        hexaFlag = True
        tetFlag = True
        fids = self.Cells[cid].Fids
        if len(fids) != 6:  # hexa need 6 faces
            hexaFlag = False
        if len(fids) != 4:
            tetFlag = False
        c_eids = set()
        c_vids = set()
        for fid in fids:
            f_eids = self.Faces[fid].Eids
            if len(f_eids) != 4:  # each face must have 4 edges
                hexaFlag = False
            if len(f_eids) != 3:  # each face must have 3 edges
                tetFlag = False
            for eid in f_eids:
                e_vids = self.Edges[eid].Vids
                c_eids.add(eid)
                for vid in e_vids:
                    if vid not in self.Faces[fid].Vids:
                        print("wrong edge vid with face ...")
                    c_vids.add(vid)
                    self.Vertices[vid].add_neighboring_vids([e_vids])
                    self.Vertices[vid].add_neighboring_eids([eid])
                    self.Vertices[vid].add_neighboring_fids([fid])
                    self.Vertices[vid].add_neighboring_cids([cid])
        self.Cells[cid].Eids = list(c_eids)
        self.Cells[cid].Vids = list(c_vids)
        if hexaFlag:
            self.Cells[cid].element_type = 6
        if tetFlag:
            self.Cells[cid].element_type = 5

    def build_element_connectivity(self):
        print("building connectivity ...")
        for vid in tqdm(self.Vertices):
            if self.Vertices[vid].removed:
                continue
            self.build_edge_connectivity(vid)
            self.build_face_connectivity(vid)
            self.build_cell_connectivity(vid)
        print("done ...")

    def build_edge_connectivity(self, vid):
        vc = self.Vertices[vid]
        v_nvids = vc.neighboring_Vids
        v_neids = vc.neighboring_Eids
        v_nfids = vc.neighboring_Fids
        v_ncids = vc.neighboring_Cids

        # all neighboring edges for a vertex are neighboring with each other
        for eid in v_neids:
            self.Edges[eid].add_neighboring_vids(v_nvids)
            self.Edges[eid].add_neighboring_eids(v_neids)

            for fid in v_nfids:
                # if edge in the face, the face is a neighboring face.
                if eid in self.Faces[fid].Eids:
                    self.Edges[eid].add_neighboring_fids([fid])

            for cid in v_ncids:
                # if edge in the cell, the cell is a neighboring face.
                if eid in self.Cells[cid].Eids:
                    self.Edges[eid].add_neighboring_cids([cid])

    def build_face_connectivity(self, vid):
        vc = self.Vertices[vid]
        v_nvids = vc.neighboring_Vids
        v_neids = vc.neighboring_Eids
        v_nfids = vc.neighboring_Fids
        v_ncids = vc.neighboring_Cids

        for fid in v_nfids:
            self.Faces[fid].add_neighboring_vids(v_nvids)
            self.Faces[fid].add_neighboring_eids(v_neids)

            f_nfids = []
            for fid_2 in v_nfids:
                if fid == fid_2:
                    continue
                eids_1 = self.Faces[fid].Eids
                eids_2 = self.Faces[fid_2].Eids
                common_edge = set.intersection(set(eids_1), set(eids_2))
                if len(common_edge) > 1:
                    self.illnessFaces.add(fid)
                    self.illnessFaces.add(fid_2)
                    # print("warning --- more than 1 common edge " + "faces : " + str(fid) + "," + str(fid_2))
                if len(common_edge) >= 1:
                    f_nfids.append(fid_2)
            self.Faces[fid].add_neighboring_fids(f_nfids)

            for cid in v_ncids:
                if fid in self.Cells[cid].Fids:
                    self.Faces[fid].add_neighboring_cids([cid])

    def build_cell_connectivity(self, vid):
        vc = self.Vertices[vid]
        v_nvids = vc.neighboring_Vids
        v_neids = vc.neighboring_Eids
        v_nfids = vc.neighboring_Fids
        v_ncids = vc.neighboring_Cids

        for cid in v_ncids:
            self.Cells[cid].add_neighboring_vids(v_nvids)
            self.Cells[cid].add_neighboring_eids(v_neids)

            c_eids = self.Cells[cid].Eids
            c_nfids = []
            for fid in v_nfids:
                if fid in self.Cells[cid].Fids:
                    continue
                f_eids = self.Faces[fid].Eids
                common_edge = set.intersection(set(f_eids), set(c_eids))
                if len(common_edge) >= 1:
                    c_nfids.append(fid)
                if len(common_edge) > 1:
                    self.illnessCells.add(cid)
                    self.illnessFaces.add(fid)
                    # print("warning --- more than 1 common edge " + "face : " + str(fid) + " - cell : " + str(cid))
            self.Cells[cid].add_neighboring_fids(c_nfids)

            c_ncids = []
            for cid_2 in v_ncids:
                if cid == cid_2:
                    continue
                fids_1 = self.Cells[cid].Fids
                fids_2 = self.Cells[cid_2].Fids
                common_face = set.intersection(set(fids_1), set(fids_2))
                if len(common_face) > 1:
                    self.illnessCells.add(cid)
                    self.illnessCells.add(cid_2)
                    # print("warning --- more than 1 common face " + "cell : " + str(cid) + "," + str(cid_2))
                if len(common_face) >= 1:
                    c_ncids.append(cid_2)
            self.Cells[cid].add_neighboring_cids(c_ncids)
            # # mesh color number equals to max number of cell neighbors
            # if len(c_ncids) + 2 > self.color_num:
            #     self.color_num = len(c_ncids) + 2

    def build(self):
        self.build_vertex_connectivity()
        self.build_element_connectivity()
        self.assign_color_to_cells()
        self.assign_color_to_edges()
        # based on face connectivity to mark is_boundary elements
        self.mark_boundary_elements()
        #
        self.mark_all_cell_type()
        #
        self.build_opposite_face_relationship()
        self.build_parallel_face_relationship()
        self.build_parallel_edge_relationship()

    def get_vertex_data(self, pid):
        return self.Vertices[pid]

    # -------------- build connectivities done -------------------
    # -------------- mark is_boundary elements ----------------------
    def mark_boundary_elements(self):
        for fid in self.Faces.keys():
            if self.Faces[fid].removed:
                continue
            neighborPolys = self.Faces[fid].neighboring_Cids
            # is_boundary face only connect to one poly
            if len(neighborPolys) != 1:
                continue
            for cid in neighborPolys:
                self.Cells[cid].is_boundary = True
            self.Faces[fid].is_boundary = True
            vids = self.Faces[fid].Vids
            for vid in vids:
                self.Vertices[vid].is_boundary = True
            eids = self.Faces[fid].Eids
            for eid in eids:
                self.Edges[eid].is_boundary = True

    # -------------- show mesh statistic -------------------
    def show_statistic(self):
        print("Mesh Path : ")
        print(self.file_path)
        print("vertex num : " + str(len(self.Vertices)))
        print("edge num : " + str(len(self.Edges)))
        print("Face num : " + str(len(self.Faces)))
        print("Cell num : " + str(len(self.Cells)))
        print("Hexa Cell num : " + str(len(self.get_hexa_cell_id_list())))
        print("Non Hexa Cell num : " + str(len(self.get_non_hexa_cell_id_list())))
        print("Hexa Ratio : " + str(len(self.get_hexa_cell_id_list()) / len(self.Cells)))
        print("Max poly Face num : " + str(self.max_poly_face_num))
        print("illness Edge num : " + str(len(self.illnessEdges)))
        print("illness Face num : " + str(len(self.illnessFaces)))
        print("illness Cell num : " + str(len(self.illnessCells)))

    def get_mesh_statistic(self):
        statList = []
        statList.append("Mesh Path " + self.file_path + "\n")
        statList.append("vertex num " + str(len(self.Vertices)) + "\n")
        statList.append("edge num " + str(len(self.Edges)) + "\n")
        statList.append("Face num " + str(len(self.Faces)) + "\n")
        statList.append("Cell num " + str(len(self.Cells)) + "\n")
        statList.append("Hexa Cell num " + str(len(self.get_hexa_cell_id_list())) + "\n")
        statList.append("Non Hexa Cell num " + str(len(self.get_non_hexa_cell_id_list())) + "\n")
        statList.append("Hexa Ratio " + str(len(self.get_hexa_cell_id_list()) / len(self.Cells)) + "\n")
        statList.append("Max poly Face num " + str(self.max_poly_face_num) + "\n")
        statList.append("illness Edge num " + str(len(self.illnessEdges)) + "\n")
        statList.append("illness Face num " + str(len(self.illnessFaces)) + "\n")
        statList.append("illness Cell num " + str(len(self.illnessCells)) + "\n")
        return statList

    def get_mesh_statistic_dict(self):
        stat_dict = {}
        stat_dict["Mesh Path"] = self.file_path
        stat_dict["vertex num"] = len(self.Vertices)
        stat_dict["edge num"] = len(self.Edges)
        stat_dict["face num"] = len(self.Faces)
        stat_dict["cell num"] = len(self.Cells)
        stat_dict["hexa cell num"] = len(self.get_hexa_cell_id_list())
        stat_dict["non hexa cell num"] = len(self.get_non_hexa_cell_id_list())
        stat_dict["hexa ratio"] = len(self.get_hexa_cell_id_list()) / len(self.Cells)
        stat_dict["max poly face num"] = self.max_poly_face_num
        return stat_dict

    # -------------- show mesh statistic done -------------------
    # -------------- opposite and dual relationship -------------------
    def build_opposite_face_relationship(self):
        # two face share an edge but do not share poly
        #  --   --
        # |f1 | f2 | 
        #  -- e1 --
        print("building opposite face relationship ...")
        for eid in self.Edges.keys():
            if self.Edges[eid].removed:
                continue
            nfs = self.Edges[eid].neighboring_Fids
            for i, fid in enumerate(nfs):
                ncs = self.Faces[fid].neighboring_Cids
                for nfid in nfs[i+1:]:
                    if fid == nfid:
                        print("still has same fid")
                        continue
                    nfncs = self.Faces[nfid].neighboring_Cids  # neighbor face's neighbor cell
                    common_cells = check_common_elements(ncs, nfncs)
                    if len(common_cells) != 0:
                        continue
                    self.Faces[fid].add_opposite_edge_face_pairs(eid, nfid)
                    self.Faces[nfid].add_opposite_edge_face_pairs(eid, fid)
        print("done ...")

    def build_parallel_face_relationship(self):
        print("building parallel face relationship ...")
        # two face in same hexa but do not share any vertex
        for cid in self.Cells.keys():
            if self.Cells[cid].removed:
                continue
            if not self.Cells[cid].is_hexa():
                continue
            fids = self.Cells[cid].Fids
            for i, fid in enumerate(fids):
                f_vids = self.Faces[fid].Vids
                for dfid in fids[i+1:]:
                    df_vids = self.Faces[dfid].Vids
                    commonVids = check_common_elements(f_vids, df_vids)
                    if len(commonVids) != 0:
                        continue
                    self.Faces[fid].add_parallel_faces(dfid)
                    self.Faces[dfid].add_parallel_faces(fid)
        print("done ...")

    def build_parallel_edge_relationship(self):
        print("building parallel edge relationship ...")
        self.clean_visited()
        # two edge in same quad face but do not share any vertex
        for eid, edge in self.Edges.items():
            if self.Edges[eid].removed:
                continue
            e_vids = edge.Vids
            edge.visited = 1
            nfs = edge.neighboring_Fids
            for fid in nfs:
                eids = self.Faces[fid].Eids
                # if face contain more than 4 edges
                # ignor it to avoid multiple dual edges
                if len(eids) != 4:
                    continue
                for deid in eids:
                    if self.Edges[deid].visited:
                        continue
                    de_vids = self.Edges[deid].Vids
                    commonVids = check_common_elements(e_vids, de_vids)
                    if len(commonVids) != 0:
                        continue
                    self.Edges[eid].add_parallels(deid)
                    self.Edges[deid].add_parallels(eid)
        self.clean_visited()
        print("done ...")

    # -------------- opposite and dual relationship done -------------------
    # -------------- return element list -------------------
    def get_vertex_list(self, skipRemoved = False):
        vertexList = []
        for vid in self.Vertices.keys():
            if skipRemoved and self.Vertices[vid].removed:
                continue
            coord = self.Vertices[vid].xyz()
            vertexList.append(coord)
        return vertexList

    def get_edge_list(self, skipRemoved = False):
        edges = []
        for eid in self.Edges.keys():
            if skipRemoved and self.Edges[eid].removed:
                continue
            edges.append(self.Edges[eid].Vids)
        return edges

    def get_face_list(self, skipRemoved = False):
        faces = []
        for fid in self.Faces.keys():
            if skipRemoved and self.Faces[fid].removed:
                continue
            faces.append(self.Faces[fid].Vids)
        return faces

    def get_face_eids_list(self, skipRemoved = False):
        faceseids = []
        for fid in self.Faces.keys():
            if skipRemoved and self.Faces[fid].removed:
                continue
            faceseids.append(self.Faces[fid].Eids)
        return faceseids

    def get_cell_list(self, skipRemoved = False):
        cells = []
        for cid in self.Cells.keys():
            if skipRemoved and self.Cells[cid].removed:
                continue
            cells.append(self.Cells[cid].Vids)
        return cells

    def get_cell_fids_list(self, skipRemoved = False):
        cellFids = []
        for cid in self.Cells.keys():
            if skipRemoved and self.Cells[cid].removed:
                continue
            fids = self.Cells[cid].Fids
            fids.sort()
            cellFids.append(fids)
        return cellFids

    def get_hybrid_cell_format_list(self, skipRemoved = False):
        # [cell fids, face orientation, cell orientation]
        cells = []
        for cid in self.Cells.keys():
            if skipRemoved and self.Cells[cid].removed:
                continue
            cells.append(self.Cells[cid].get_hybrid_cell_format())
        return cells

    def get_hexa_cell_id_list(self, skipRemoved = False):
        hexaCellIds = []
        for cid in self.Cells.keys():
            if skipRemoved and self.Cells[cid].removed:
                continue
            if self.Cells[cid].is_hexa():
                hexaCellIds.append(cid)
        return hexaCellIds

    def get_non_hexa_cell_id_list(self, skipRemoved = False):
        self.mark_all_cell_type()
        nonHexaCellIds = []
        for cid in self.Cells.keys():
            if skipRemoved and self.Cells[cid].removed:
                continue
            if not self.Cells[cid].is_hexa():
                nonHexaCellIds.append(cid)
        return nonHexaCellIds

    def get_boundary_vids(self, skipRemoved = False):
        boundaryVids = set()
        for fid in self.Faces.keys():
            if skipRemoved and self.Faces[fid].removed:
                continue
            ncells = self.Faces[fid].neighboring_Cids
            if len(ncells) == 1:
                boundaryVids.update(self.Faces[fid].Vids)
        return list(boundaryVids)

    def get_boundary_eids(self, skipRemoved = False):
        boundary = set()
        for eid in self.Edges:
            edge = self.Edges[eid]
            if edge.is_boundary:
                boundary.add(eid)
        return boundary

    def get_boundary_fids(self, skipRemoved = False):
        boundary = set()
        for fid in self.Faces:
            face = self.Faces[fid]
            if face.is_boundary:
                boundary.add(fid)
        return boundary

    def get_edge_neighbor_bnd_fids(self, eid):
        edge = self.Edges[eid]
        bnd_faces = set()
        for nfid in edge.neighboring_Fids:
            face = self.Faces[nfid]
            if face.is_boundary:
                bnd_faces.add(nfid)
        return bnd_faces

    def remove_not_closed_boundary_fids(self, boundary_fids):
        # boundary_fids = self.get_boundary_fids()
        face_stack = list(boundary_fids)
        removed_fids = set()
        while face_stack:
            fid = face_stack.pop()
            if fid in removed_fids:
                continue
            face = self.Faces[fid]
            remove_face = False
            for eid in face.Eids:
                bnd_fids = self.get_edge_neighbor_bnd_fids(eid)
                bnd_fids = set(bnd_fids) - removed_fids
                if len(bnd_fids) < 2:
                    removed_fids.add(fid)
                    remove_face = True
                    continue
            if not remove_face:
                continue
            face_stack.extend(face.neighboring_Fids)
        real_bnd_fids = set()
        for fid in boundary_fids:
            if fid in removed_fids:
                continue
            real_bnd_fids.add(fid)
        return real_bnd_fids, removed_fids

    def get_invalid_bnd_fids(self, bnd_fid):
        face = self.Faces[bnd_fid]
        removed_fids = set()
        if not face.is_boundary:
            return removed_fids
        if len(face.Eids) > 4:
            return removed_fids
        for nfid in face.neighboring_Fids:
            n_face = self.Faces[nfid]
            if not n_face.is_boundary:
                continue
            if len(n_face.Eids) > 4:
                continue
            shared_eids = check_common_elements(n_face.Eids, face.Eids)
            if len(shared_eids) > 1:
                removed_fids.add(nfid)
        if removed_fids:
            removed_fids.add(bnd_fid)
        return removed_fids

    # two bnd face can not share more than one bnd edge
    # only for quad and tri faces
    def remove_fold_bnd_fids(self, bnd_fids):
        removed_fids = set()
        face_stack = list(bnd_fids)
        while face_stack:
            fid = face_stack.pop()
            if fid in removed_fids:
                continue
            invalid_fids = self.get_invalid_bnd_fids(fid)
            removed_fids.update(invalid_fids)
        real_bnd_fids = set()
        for fid in bnd_fids:
            if fid in removed_fids:
                continue
            real_bnd_fids.add(fid)
        return real_bnd_fids, removed_fids

    def get_manifold_boundary_fids(self):
        boundary_fids = self.get_boundary_fids()
        real_bnd_fids = set()
        not_cleaned_bnd = True
        n = 0
        while not_cleaned_bnd:
            print(f"cleaning bnd faces ... iter {n}", end="\r")
            removed_fids = set()
            boundary_fids, removed_unfold_bnds = self.remove_fold_bnd_fids(boundary_fids)
            removed_fids.update(removed_unfold_bnds)
            boundary_fids, removed_not_closed_bnds = self.remove_not_closed_boundary_fids(boundary_fids)
            removed_fids.update(removed_not_closed_bnds)
            if len(removed_fids) == 0:
                not_cleaned_bnd = False
            n += 1
        return boundary_fids

    # -------------- return element list done -------------------
    # -------------- check manifold conditions -------------------
    def check_face_manifold(self):
        # non manifold face:
        #   1. two neighbor faces share more than one edge
        #   2. a face included in more than two cells
        #   3. two face have more than two common vids but not neighboring with each other
        #       3.1 hard to check manifold for #3 if model has a hole inside volumn
        #       3.2 we need to assume no gap inside volumn
        print("checking face manifold ...")
        nonManifoldFids = []
        for fid in self.Faces.keys():
            manifold = True
            face = self.Faces[fid]
            f_vids = face.Vids
            f_eids = face.Eids
            nfids = face.neighboring_Fids
            # 1. two neighbor faces share more than on edge
            for nfid in nfids:
                nface = self.Faces[nfid]
                nf_eids = nface.Eids
                commonEdges = check_common_elements(f_eids, nf_eids)
                if len(commonEdges) != 1:
                    manifold = False
                    print("face : " + str(fid) + " and its neighbor face : " + str(nfid) + " have " + str(len(commonEdges)) + " common edges.")
            # 2. face included in more than two cells, or 0 cell
            ncids = len(set(face.neighboring_Cids))
            if ncids == 0 or ncids > 2:
                manifold = False
                print("face : " + str(fid) + " in " + str(ncids) + " cells.")
            # 3. two face have more than one common vids but not neighboring with each other
            # for otherFid in self.Faces.keys():
            #     if otherFid == fid or otherFid in nfids:
            #         continue
            #     otherFace = self.Faces[otherFid]
            #     of_vids = otherFace.Vids
            #     commonVids = check_common_elements(f_vids, of_vids)
            #     if len(commonVids) > 1:
            #         manifold = False
            #         print("face : " + str(fid) + " and face : " + str(otherFid) + " have " + str(len(commonVids)) + " common vids.")
            if not manifold:
                nonManifoldFids.append(fid)
        print("done ...")
        return nonManifoldFids

    def check_edge_manifold(self):
        # non manifold edge : 
        #   1. a edge only be included in one face
        #   2. two edges has same vids
        #   3. a edge has two same vids
        print("checking edge manifold ...")
        nonManifoldEids = []
        for eid in self.Edges.keys():
            manifold = True
            edge = self.Edges[eid]
            if len(set(edge.Vids)) != 2:
                manifold = False
                print("edge : " + str(eid) + " does not have two unique vids.")
            if len(set(edge.neighboring_Fids)) == 1:
                manifold = False
                print("edge : " + str(eid) + " only in one face.")
            neids = edge.neighboring_Eids
            for neid in neids:
                nedge = self.Edges[neid]
                if set(edge.Vids) == set(nedge.Vids):
                    manifold = False
                    print("edge : " + str(eid) + " and edge : " + str(neid) + " have same vids.")
            if not manifold:
                nonManifoldEids.append(eid)
        print("done ...")
        return nonManifoldEids

    def check_cell_manifold(self):
        # non manifold cells:
        #   1. two cells contain more than one common face
        print("checking cell manifold ...")
        nonManifoldCids = []
        for cid in self.Cells.keys():
            manifold = True
            cell = self.Cells[cid]
            fids = cell.Fids
            for otherCid in self.Cells.keys():
                if cid == otherCid:
                    continue
                otherCell = self.Cells[otherCid]
                oc_fids = otherCell.Fids
                commonFids = check_common_elements(fids, oc_fids)
                if len(commonFids) > 1:
                    manifold = False
                    print("cell : " + str(cid) + " and cell : " + str(otherCid) + " have more than two common face.")
            if not manifold:
                nonManifoldCids.append(cid)
        print("done ...")
        return nonManifoldCids

    def check_vertex_manifold(self):
        # non manifold vertex:
        #   1. a vertex does not have two edges in a face
        #   2. vertex not included in any face
        print("checking vertex manifold ...")
        nonManifoldVidPairs = []
        for vid in self.Vertices.keys():
            vertex = self.Vertices[vid]
            neighborFids  = vertex.neighboring_Fids
            manifold = True
            for fid in neighborFids:
                neighborV = self.get_all_neighbor_vids_for_vid_in_face(fid, vid)
                if len(neighborV) != 2:
                    print("this vertex : " + str(vid) + " contains {} edges in the face : ".format(len(neighborV)) + str(fid))
                    print("face vids : ")
                    print(self.Faces[fid].Vids)
                    print("neighborV : ")
                    print(neighborV)
                    print("neighbor face vids : ")
                    for fff in neighborFids:
                        print(self.Faces[fff].Vids)
                    # currently ignor the invalid face
                    for v in neighborV:
                        pair = [vid, v]
                        pair.sort()
                        # slow insert
                        if pair not in nonManifoldVidPairs:
                            nonManifoldVidPairs.append(pair)
        print("done ...")
        return list(nonManifoldVidPairs)

    # -------------- check manifold conditions done -------------------
    # -------------- find edges contain vid in a face -------------------
    def get_all_neighbor_vids_for_vid_in_face(self, fid, vid):
        ans = []
        face = self.Faces[fid]
        vertex = self.Vertices[vid]
        neighborV = vertex.neighboring_Vids
        vids = face.Vids
        for nvid in vids:
            if nvid in neighborV:
                ans.append(nvid)
        return ans

    # -------------- find edges contain vid in a face done -------------------
    def is_face_on_domain_boundary(self, cids, fid):
        ncids = check_common_elements(self.Faces[fid].neighboring_Cids, cids)
        if len(ncids) > 1:
            return False
        return True

    def cids_to_fids(self, cids, only_boundary=False):
        all_fids = set()
        for cid in cids:
            for fid in self.Cells[cid].Fids:
                if only_boundary and not self.is_face_on_domain_boundary(cids, fid):
                    continue
                all_fids.add(fid)
        return list(all_fids)

    # -------------- give cell id return face vids -------------------
    def cids_to_face_vids(self, cids, only_boundary=False):
        faceVids = []
        for cid in cids:
            fids = self.Cells[cid].Fids
            for fid in fids:
                if only_boundary and not self.is_face_on_domain_boundary(cids, fid):
                    continue
                vids = self.Faces[fid].Vids
                faceVids.append(vids)
        return faceVids

    def fids_to_face_vids(self, fids):
        faceVids = []
        for fid in fids:
            vids = self.Faces[fid].Vids
            faceVids.append(vids)
        return faceVids

    def eids_to_edge_vids(self, eids):
        edgeVids = []
        for eid in eids:
            vids = self.Edges[eid].Vids
            edgeVids.append(vids)
        return edgeVids

    # -------------- give cell id return face vids done -------------------
    def find_parallel_fid_in_cell(self, cid, fid):
        hexa = self.Cells[cid]
        hfids = hexa.Fids
        parallel_fid = -1 # the fid that parallel with given fid
        for pfid in self.Faces[fid].parallels:
            if pfid in hfids:
                parallel_fid = pfid
        # if could not find parallel fid
        # return the given -1
        return parallel_fid

    def edge_vids_between_two_fids_in_cell(self, fid_0, fid_1):
        face_0 = self.Faces[fid_0]
        face_1 = self.Faces[fid_1]
        allVids = face_0.Vids + face_1.Vids
        edgeVidMap = {vid:[] for vid in allVids}
        for neid in face_0.neighboring_Eids:
            if neid in face_1.Eids:
                continue
            vids = self.Edges[neid].Vids
            for vid in vids:
                if vid in face_1.Vids:
                    edgeVidMap[vids[0]].append(vids[1])
                    edgeVidMap[vids[1]].append(vids[0])
        return edgeVidMap

    # -------------- find element index based on vids -------------------
    def get_eid_based_vids(self, evids):
        if hasattr(self, "edgeVids"):
            if len(self.edgeVids) != len(self.Edges):
                self.edgeVids = self.get_edge_list()
        else:
            self.edgeVids = self.get_edge_list()
        evids.sort()
        try:
            eid = self.edgeVids.index(evids)
        except:
            eid = -1
        return eid

    def get_fid_based_vids(self, fvids):
        # face vids is not sorted.
        if hasattr(self, "faceVids"):
            if len(self.faceVids) != len(self.Faces):
                self.faceVids = self.get_face_list()
        else:
            self.faceVids = self.get_face_list()
        try:
            fid = self.faceVids.index(fvids)
        except:
            fid = -1
        return fid

    def get_cid_based_vids(self, cvids):
        # face vids is not sorted.
        if hasattr(self, "cellVids"):
            if len(self.cellVids) != len(self.Cells):
                self.cellVids = self.get_cell_list()
        else:
            self.cellVids = self.get_cell_list()
        cvids.sort() # maybe necessary
        try:
            cid = self.cellVids.index(cvids)
        except:
            cid = -1
        return cid

    def get_cid_based_fids(self, cfids):
        # face vids is not sorted.
        if hasattr(self, "cellFids"):
            if len(self.cellFids) != len(self.Cells):
                self.cellFids = self.get_cell_fids_list()
        else:
            self.cellFids = self.get_cell_fids_list()
        cfids.sort() # maybe necessary
        try:
            cid = self.cellFids.index(cfids)
        except:
            cid = -1
        return cid

    # -------------- find element index based on vids done -------------------
    # -------------- update element attribute -------------------
    def update_face_attribute(self, fid):
        # update a exiting face's edge attribute
        face = self.Faces[fid]
        fvids = face.Vids
        edgesVids, _ = face_list_to_edge_list([fvids])
        newEids = set()
        for evids in edgesVids:
            eid = self.get_eid_based_vids(evids)
            newEids.add(eid)
        face.update_eids(list(newEids))

    def update_cell_attribute(self, cid):
        # cell must contain fids
        # assume faces are correct
        cell = self.Cells[cid]
        fids = cell.Fids
        # new eids
        newEids = set()
        # new vids # not ordered
        newVids = set()
        for fid in fids:
            face = self.Faces[fid]
            eids = face.Eids
            newEids.update(eids)
            vids = face.Vids
            newVids.update(vids)
        cell.Eids = list(newEids)
        cell.Vids = list(newVids)

    # -------------- update all element attributes -------------------
    # face and cell element attributes are needed to be updated
    def update_faces_attribute(self, fidMap):
        faceMapValues = fidMap.values()
        for fid in set(faceMapValues):
            if fid < 0:
                continue
            self.update_face_attribute(fid)

    def update_cells_attribute(self, cidMap):
        cellMapValues = cidMap.values()
        for cid in set(cellMapValues):
            if cid < 0:
                continue
            self.update_cell_attribute(cid)

    # -------------- update all element attributes done -------------------
    def get_vert_neighbor_cell_groups(self, vid):
        vert = self.Vertices[vid]
        v_ncids = vert.neighboring_Cids
        g = []
        for cid in v_ncids:
            cell = self.Cells[cid]
            in_scope_ncids = check_common_elements(cell.neighboring_Cids, v_ncids)
            in_scope_ncids.add(cid)
            g.append(in_scope_ncids)
        g = reduce_id_groups(g)
        return g

    def mark_non_manifold_vertex(self):
        print("detecting non manifold vertices ...")
        for vid, vert in self.Vertices.items():
            # vert cell group
            g = self.get_vert_neighbor_cell_groups(vid)
            if len(g) != 1:
                vert.is_non_manifold = True
            else:
                vert.is_non_manifold = False
            if vert.is_non_manifold:
                continue
            # vert in cell face number
            vert_nfids = vert.neighboring_Fids
            for cid in vert.neighboring_Cids:
                cell = self.Cells[cid]
                if len(check_common_elements(cell.Fids, vert_nfids)) < 3:
                    vert.is_non_manifold = True
                    break
                else:
                    vert.is_non_manifold = False
        print("detecting non manifold vertices ... done")

    # -------------- mark singularity elements ------------------------
    # all non hexa cell's edges are also marked as singular edges
    def mark_singularity_elements(self):
        self.mark_all_cell_type()
        self.mark_non_manifold_vertex()
        for eid, edge in self.Edges.items():
            nc = len(edge.neighboring_Cids)
            if edge.is_boundary:
                if nc != Boundary_Regular_E:
                    edge.is_singular = True
                else:
                    edge.is_singular = False
            elif nc != Interior_Regular_E:
                edge.is_singular = True
            else:
                edge.is_singular = False
            # non hexa cell's edge are also singular
            for c in edge.neighboring_Cids:
                if self.Cells[c].is_hexa():
                    continue
                edge.is_singular = True
                # count non hexa edge and hexa singular edge separately.
                self.non_hexa_eids.add(eid)
            # count non hexa edge and hexa singular edge separately.
            if (eid not in self.non_hexa_eids) and edge.is_singular:
                self.hexa_singular_eids.add(eid)
            # non manifold vert - neighboring edges
            if not edge.is_singular:
                for vid in edge.Vids:
                    if self.Vertices[vid].is_non_manifold:
                        edge.is_singular = True
                        break
            # mark all vertices of singular edge as singular vertex
            if edge.is_singular:
                for evid in edge.Vids:
                    self.Vertices[evid].is_singular = True

    # -------------- mark singularity elements done -------------------
    # -------------- find a hexa box that surround non hexa elements -------------------
    def get_cids_for_hexa_box(self):
        row_hexa_box = set()
        for eid, edge in self.Edges.items():
            all_hexa_flag = True
            for ncid in edge.neighboring_Cids:
                cell = self.Cells[ncid]
                if cell.is_hexa():
                    continue
                all_hexa_flag = False
            if not all_hexa_flag:
                row_hexa_box.update(edge.neighboring_Cids)
        return row_hexa_box

    def get_singular_edges_with_pure_hexa_neighbor(self):
        all_hexa_neighbor_singular_edge = set()
        for eid, edge in self.Edges.items():
            all_hexa = True
            for cid in edge.neighboring_Cids:
                cell = self.Cells[cid]
                if cell.is_hexa():
                    continue
                all_hexa = False
            if not all_hexa:
                continue
            valence = len(edge.neighboring_Cids)
            if edge.is_boundary:
                if valence != 2:
                    all_hexa_neighbor_singular_edge.add(eid)
            else:
                if valence != 4:
                    all_hexa_neighbor_singular_edge.add(eid)
        return all_hexa_neighbor_singular_edge

    # -------------- find a hexa box that surround non hexa elements done -------------------
    def find_boundary_face_of_non_hexa_cells(self):
        non_hexa_cell_boundary_faces = set()
        for cid, cell in self.Cells.items():
            if cell.is_hexa():
                continue
            for f in cell.Fids:
                face = self.Faces[f]
                if face.is_boundary:
                    non_hexa_cell_boundary_faces.add(f)
        return non_hexa_cell_boundary_faces

    def get_isolated_vert_in_non_hexa_cells(self):
        isolated_vert = set()
        for cid, cell in self.Cells.items():
            for vid in cell.Vids:
                vert = self.Vertices[vid]
                common_eids = check_common_elements(cell.Eids, vert.neighboring_Eids)
                if len(common_eids) == len(vert.neighboring_Eids):
                    isolated_vert.add(vid)
        return isolated_vert

    def get_non_hexa_cids_that_neighbor_hexa(self):
        non_hexa_cids = set()
        for cid, cell in self.Cells.items():
            if cell.is_hexa():
                continue
            for ncid in cell.neighboring_Cids:
                if self.Cells[ncid].is_hexa():
                    non_hexa_cids.add(cid)
                    break
        return non_hexa_cids

    # edges surround by non-hex
    def get_local_eids(self):
        local_eids = []
        for eid, edge in self.Edges.items():
            local = True
            for ncid in edge.neighboring_Cids:
                cell = self.Cells[ncid]
                if cell.is_hexa():
                    local = False
                    break
            if local:
                local_eids.append(eid)
        return local_eids

    # if a face contain a neighboring face that in more than one cids
    # highligt those faces
    def local_fids(self):
        local_fids = []
        for fid, face in self.Faces.items():
            local = False
            for nfid in face.neighboring_Fids:
                nface = self.Faces[nfid]
                common_cids = check_common_elements(face.neighboring_Cids, nface.neighboring_Cids)
                if len(common_cids) > 1:
                    local = True
            if local:
                local_fids.append(fid)
        return local_fids

    def get_hanging_fids(self):
        hanging_fids = []
        for fid, face in self.Faces.items():
            for eid in face.Eids:
                edge = self.Edges[eid]
                if len(edge.neighboring_Fids) <= 1:
                    if fid in hanging_fids:
                        continue
                    hanging_fids.append(fid)
        return hanging_fids

    # -------------------------------------------------------------
    # opposite element extraction
    # -------------------------------------------------------------
    def trace_patch_faces(self, startFid):
        fidStack = [startFid]
        patchFids = set(fidStack)
        while len(fidStack) > 0:
            cur_fid = fidStack.pop()
            if self.Faces[cur_fid].visited:
                continue
            self.Faces[cur_fid].visited = True
            op_fids = self.get_opposite_faces_of_face(cur_fid)
            if len(op_fids) > 0:
                fidStack.extend(op_fids)
                patchFids.update(op_fids)
        return patchFids

    def get_opposite_faces_of_face(self, fid):
        oppositeFaces = []
        face = self.Faces[fid]
        for eid in face.Eids:
            oppositeFid = self.get_opposite_face_on_regular_eid(fid, eid)
            if len(oppositeFid) > 0:
                oppositeFaces.extend(oppositeFid)
        return oppositeFaces

    def get_opposite_face_on_regular_eid(self, fid, eid):
        # singularity edge should only contain one opposite face of an edge at a time
        oppositeFace = []
        edge = self.Edges[eid]
        # stop at the singular edge
        # if edge.is_singular or edge.is_boundary:
        #     return oppositeFace
        if edge.is_singular or edge.is_boundary:
            return oppositeFace
        oppositeFace = self.get_opposite_face_on_any_type_eid(fid, eid)
        # only contain one opposite face
        return oppositeFace

    def get_opposite_face_on_any_type_eid(self, fid, eid):
        # singularity edge should only contain one opposite face of an edge at a time
        oppositeFace = []
        face = self.Faces[fid]
        edge = self.Edges[eid]
        for nfid in edge.neighboring_Fids:
            nface = self.Faces[nfid]
            if nface.visited:
                continue
            # ToDo : require two faces are same type
            # if face.is_boundary != nface.is_boundary:
            #     continue
            sharingCids = check_common_elements(face.neighboring_Cids, nface.neighboring_Cids)
            # if both face and nface contain same neighboring cids
            # they are at the same cid
            if len(sharingCids) > 0:
                continue
            oppositeFace.append(nfid)
            # if a regular edge only has one opposite face
            # can return here
            # return nfid
        if len(oppositeFace) > 1:
            # print(edge.neighboring_Cids)
            self.multi_opposite_fids.add(fid)
            self.multi_opposite_fids.update(oppositeFace)
            # print(oppositeFace)
            print("a regular edge contain a face have two opposite face....")
            # return []
        # only contain one opposite face
        return oppositeFace

    # -------------- is two edge in same face -------------------
    def is_two_edge_in_same_face(self, eid_0, eid_1):
        edge_0 = self.Edges[eid_0]
        edge_1 = self.Edges[eid_1]
        nf_0 = edge_0.neighboring_Fids
        nf_1 = edge_1.neighboring_Fids
        commonFids = check_common_elements(nf_0, nf_1)
        if len(commonFids) > 0:
            return True
        return False

    def is_two_edge_in_same_cell(self, eid_0, eid_1):
        edge_0 = self.Edges[eid_0]
        edge_1 = self.Edges[eid_1]
        nc_0 = edge_0.neighboring_Cids
        nc_1 = edge_1.neighboring_Cids
        commonCids = check_common_elements(nc_0, nc_1)
        if len(commonCids) > 0:
            return True
        return False

    # https://github.com/Cotrik/CotrikMesh/blob/d890c9e596fa3242342b8a1b7fc7b9ce194f03a7/libcotrik/src/BaseComplex.cpp#L1867
    def get_next_edge(self, vid, eid):
        nextEids = set()
        vert = self.Vertices[vid]
        neids = vert.neighboring_Eids
        for neid in neids:
            if self.is_two_edge_in_same_face(eid, neid):
                continue
            # if self.is_two_edge_in_same_cell(eid, neid):
            #     continue
            # ToDo : consider if a edge trace from bnd to internal
            # ToDo : or trace from internal to bnd
            # if not self.edge_end_vertices_are_same_type(neid):
            #     continue

            nextEids.add(neid)
        # if len(nextEids) > 1:
        #     print(f'edge {eid} contain {len(nextEids)} edges connected by vert {vid}')
        # if len(nextEids) == 0:
        #     print("not next eids is found")
        return list(nextEids)
    # -------------- is two edge in same face done -------------------

    # next edge require next edge contain same type
    # both internal or both boundary
    def edge_end_vertices_are_same_type(self, eid):
        same_flag = False
        target_flag = False
        if len(self.Edges[eid].Vids) != 2:
            return False
        for vid in self.Edges[eid].Vids:
            vert = self.Vertices[vid]
            if vert.is_boundary:
                same_flag = not same_flag
        if target_flag == same_flag:
            return True
        else:
            return False

    def is_two_fids_non_comformal(self, fid_0, fid_1):
        non_comformal = False
        if fid_1 == fid_0:
            return non_comformal
        face_0 = self.Faces[fid_0]
        face_1 = self.Faces[fid_1]
        common_vids = check_common_elements(face_0.Vids, face_1.Vids)
        if len(common_vids) < 2:
            return non_comformal
        # check whether face 0 contain a eid make non comformal neighbor
        for eid in face_0.Eids:
            if eid in face_1.Eids:
                continue
            edge = self.Edges[eid]
            in_face_vids = check_common_elements(edge.Vids, face_1.Vids)
            if len(in_face_vids) == len(edge.Vids):
                non_comformal = True
        # check whether face 0 contain a eid make non comformal neighbor
        for eid in face_1.Eids:
            if eid in face_0.Eids:
                continue
            edge = self.Edges[eid]
            in_face_vids = check_common_elements(edge.Vids, face_0.Vids)
            if len(in_face_vids) == len(edge.Vids):
                non_comformal = True
        return non_comformal

    def get_non_comformal_fids(self):
        non_comformal_fids = set()
        for vid, vert in self.Vertices.items():
            neighbor_fids = vert.neighboring_Fids
            for nfid in neighbor_fids:
                for compare_fid in neighbor_fids:
                    if self.is_two_fids_non_comformal(nfid, compare_fid):
                        non_comformal_fids.update([nfid, compare_fid])
        return list(non_comformal_fids)

    def get_ruler_characteristic(self, cid):
        cell = self.Cells[cid]
        return len(cell.Fids) - len(cell.Eids) + len(cell.Vids)

    def is_cell_homeomorphism_to_sphere(self, cid):
        ruler = self.get_ruler_characteristic(cid)
        if ruler != 2:
            return False
        else:
            return True

    def identify_cells_homeomorphism(self):
        homeomorphism_to_sphere = set()
        not_homeomorphism = set()
        for cid in self.Cells:
            if self.get_ruler_characteristic(cid):
                homeomorphism_to_sphere.add(cid)
            else:
                not_homeomorphism.add(cid)
        return homeomorphism_to_sphere, not_homeomorphism

    def is_hexa(self, cid):
        cell = self.Cells[cid]
        hexa = True
        if len(cell.Fids) != 6:
            hexa = False
        elif len(cell.Vids) != 8:
            hexa = False
        elif len(cell.Eids) != 12:
            hexa = False
        else:
            for fid in cell.Fids:
                cf = self.Faces[fid]
                if not cf.is_quad():
                    hexa = False
                    break
        if hexa:
            cell.element_type = 6
            cell.hex_flag = True
        else:
            cell.element_type = 4

        return hexa

    def get_cell_parallel_groups(self, cid):
        cell = self.Cells[cid]
        parallel_groups = []
        for eid in cell.Eids:
            edge = self.Edges[eid]
            in_cell_parallel_edges = check_common_elements(edge.parallels, cell.Eids)
            # print(f"------ in cell : before {in_cell_parallel_edges}")
            in_cell_parallel_edges.add(eid)
            # print(f"------ in cell : after {in_cell_parallel_edges}")
            parallel_groups.append(in_cell_parallel_edges)
        parallel_groups = reduce_id_groups(parallel_groups)
        return parallel_groups

    def mark_all_cell_type(self):
        # cell type of hexa is 6 (hexa)
        # other type of cell is 4 (cell)
        for cid in self.Cells:
            res = self.is_hexa(cid)
            if res:
                continue
            # update face hex_flag
            for fid in self.Cells[cid].Fids:
                if self.Faces[fid].hex_flag:
                    self.Faces[fid].hex_flag = res
            # update edge hex_flag
            for eid in self.Cells[cid].Eids:
                if self.Edges[eid].hex_flag:
                    self.Edges[eid].hex_flag = res
            # update vertex hex_flag
            for vid in self.Cells[cid].Vids:
                if self.Vertices[vid].hex_flag:
                    self.Vertices[vid].hex_flag = res

    def get_edge_length(self, eid):
        vids = self.Edges[eid].Vids
        v_0_xyz = np.array(self.Vertices[vids[0]].xyz())
        v_1_xyz = np.array(self.Vertices[vids[1]].xyz())
        distance = np.linalg.norm(v_0_xyz - v_1_xyz)
        return distance

    def get_average_edge_length(self):
        sum_edge_length = 0
        for eid in self.Edges:
            sum_edge_length += self.get_edge_length(eid)
        return sum_edge_length / len(self.Edges)

    def assign_color_to_cells(self):
        for cid in self.Cells:
            neighbor_colors = set()
            for ncid in self.Cells[cid].neighboring_Cids:
                cc = self.Cells[ncid]
                neighbor_colors.add(cc.color)
            if len(neighbor_colors) >= self.color_num - 2:
                self.color_num = len(neighbor_colors) + 2  # extra one for the input cell
            candidate_colors = set()
            for i in range(self.color_num):
                # keep rad color reserved
                if i == 0:
                    continue
                if i not in neighbor_colors:
                    candidate_colors.add(i)
            picked_color = random.choice(list(candidate_colors))
            # pickColor = list(candidate_colors)[int(len(candidate_colors) / 2)]
            self.Cells[cid].color = picked_color

    def assign_color_to_edges(self):
        # determine color number
        for eid in self.Edges:
            neids = len(self.Edges[eid].neighboring_Eids)
            if neids + 1 > self.edge_color_num:
                self.edge_color_num = neids + 1
        # assgin color
        for eid, edge in self.Edges.items():
            neighbor_colors = set()
            for neid in edge.neighboring_Eids:
                ne = self.Edges[neid]
                neighbor_colors.add(ne.color)
            candidate_colors = set()
            for i in range(self.edge_color_num):
                # keep rad color reserved
                if i == 0:
                    continue
                if i not in neighbor_colors:
                    candidate_colors.add(i)
            picked_color = random.choice(list(candidate_colors))
            edge.color = picked_color

    # hex cells in a groups are connected
    def get_hex_cell_groups_around_edge(self, eid):
        edge = self.Edges[eid]
        cells_groups = []
        for cid in edge.neighboring_Cids:
            cell = self.Cells[cid]
            if not cell.is_hexa():
                continue
            tmp_cids = check_common_elements(cell.neighboring_Cids, edge.neighboring_Cids)
            ncids = set()
            for tc in tmp_cids:
                tmp = self.Cells[tc]
                if tmp.is_hexa():
                    ncids.add(tc)
            ncids.add(cid)
            cells_groups.append(ncids)
        cells_groups = reduce_id_groups(cells_groups)
        return cells_groups

    def check_duplicate_faces(self):
        existing_faces = []
        duplicate_count = 0
        for fid in self.Faces:
            face = self.Faces[fid]
            vids = set(face.Vids)
            if vids in existing_faces:
                duplicate_count += 1
                continue
            existing_faces.append(vids)
        if duplicate_count != 0:
            print(f"mesh contain {duplicate_count} duplicate faces.")
        return existing_faces

    def check_duplicate_vertices(self):
        xyz_list = self.get_vertex_list()
        duplicate_count = find_duplicate_elements(xyz_list)
        if duplicate_count != 0:
            print(f"mesh contain {duplicate_count} duplicate vertices.")
        return
