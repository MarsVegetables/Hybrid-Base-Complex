'''
Author: Lei Si
Date: 2022-08-15 20:46:25
LastEditTime: 2022-08-24 16:58:56
LastEditors: Lei Si
Description: 
FilePath: /hexa_dominant_layer_simplification/meshLoader/hybrid_file_loader.py
YooooHooooo
'''
import copy

from mesh_structure.multi_level_hybrid_mesh import MultiLevelHybridMesh
from mesh_structure.basic_mesh import BasicMesh
from mesh_structure.mesh_common import add_face_cell_orientation, check_common_elements, varify_and_clean_file_data
import pathlib


class MeshLoader:
    def __init__(self, file_path, load_as_hybrid=False):
        self.file_name = pathlib.Path(file_path).stem
        self.file_ext = pathlib.Path(file_path).suffix
        self.input_file_path = file_path
        # True : when the structure includes component and supercell
        # False: when the structure is used only for load a file
        self.hybrid_structure = load_as_hybrid
        self.mesh = []

    def load(self):
        vertices = []
        faces = []
        polys = []
        if self.hybrid_structure:
            self.mesh = MultiLevelHybridMesh(self.input_file_path, "Hybrid")
        else:
            self.mesh = BasicMesh(self.input_file_path, "Basic")

        if self.file_ext == ".mesh":
            vertices, faces, polys = load_mesh_mesh(self.input_file_path)
        elif self.file_ext == ".HYBRID":
            vertices, faces, polys = load_mesh_hybrid(self.input_file_path)
            faces, polys = varify_and_clean_file_data(faces, polys)
        elif self.file_ext == ".hedra":
            vertices, faces, polys = load_mesh_hedra(self.input_file_path)
            # faces, polys = varify_and_clean_file_data(faces, polys)

        if vertices and faces and polys:
            self.mesh.build_mesh_from_vfc(vertices, faces, polys)
            return self.mesh
        else:
            return []


# input : cell list
# cell list structure
#   [[face id list],
#    [face orientation list],
#    cell orientation]
# if value of orientation is 1, the face or cell need to change the normal direction.
def load_mesh_hybrid(f_path):
    with open(f_path, 'r') as my_file:
        nv, nf, nc = [int(i) for i in my_file.readline().split()]
        vertices = [[float(i) for i in my_file.readline().split()] for j in range(nv)]
        faces = [[int(i) for i in my_file.readline().split()] for j in range(nf)]
        # nc number include poly face list, poly face orientation list, poly orientation list
        poly_num = nc / 3
        poly = [[[int(i) for i in my_file.readline().split()], [int(j) for j in my_file.readline().split()], 0] for k in range(int(poly_num))]
        for k in range(int(poly_num)):
            poly[k][2] = int(my_file.readline())
        return vertices, faces, poly


# hydra file format is defined in
# https://github.com/mlivesu/cinolib/blob/master/include/cinolib/io/read_HEDRA.cpp
def load_mesh_mesh(f_path):
    # load a hybrid mesh from .mesh file
    with open(f_path, 'r') as my_file:
        vertices, faces, newPoly = [], [], []
        formatVersion = my_file.readline()
        dim = my_file.readline()
        element = str(my_file.readline())
        if element == "Vertices\n":
            nv = int(my_file.readline())
        else:
            print("first element is not vertices ...")
            return vertices, faces, newPoly
        # skip last 0 value
        vertices = [[float(i) for i in my_file.readline().split()[:-1]] for j in range(nv)]
        element = str(my_file.readline())
        if element == "Hexahedra\n":
            poly_num = int(my_file.readline())
        else:
            print("second element is not hexa ...")
            return vertices, faces, newPoly
        poly = [[int(i) - 1 for i in my_file.readline().split()[:-1]] for j in range(poly_num)]
        faces, polyFids = covert_mesh_poly_to_face(poly)
        newPoly = add_face_cell_orientation(polyFids)
        return vertices, faces, newPoly


def load_mesh_hedra(f_path):
    # load a hybrid mesh from .mesh file
    with open(f_path, 'r') as my_file:
        nv, nf, nc = [int(i) for i in my_file.readline().split()]
        vertices = [[float(i) for i in my_file.readline().split()] for j in range(nv)]
        faces = []
        for j in range(nf):
            cur_line = my_file.readline().split()
            nvids = cur_line[0]
            # cur_face = [int(nvids)]
            cur_face = []
            for v in cur_line[1:]:
                cur_face.append(int(v) - 1)
            faces.append(cur_face)
        polyFids = []
        for k in range(nc):
            cur_line = my_file.readline().split()
            nfids = cur_line[0]
            cur_fids = []
            for f in cur_line[1:]:
                cur_fids.append(abs(int(f)) - 1)
            polyFids.append(cur_fids)
        newPoly = add_face_cell_orientation(polyFids)
        return vertices, faces, newPoly


# /*
#  * HEXALAB vertex and face orderning
#  *
#  *      6-------7       +   -   +       +-------+       +-------+
#  *     /|      /|      /|      /|      /   3   /        |   5   |
#  *    2-------3 |     + | -   + |     +-------+       +-------+ |
#  *    | |     | |     |0|     |1|                     | |     | |
#  *    | 4-----|-5     | +   - | +       +-------+     | +-----|-+
#  *    |/      |/      |/      |/       /   2   /      |   4   |
#  *    0-------1       +   -   +       +-------+       +-------+
#  *
#  *      Vertex                           Faces
#  *
#  *
#  * 2---6---7                               +-0-+-2-+
#  * | 0 | 5 |                               3   1   3
#  * 0---4---5---7         0---1  +-0-+      +-2-+-0-+-2-+
#  *     | 2 | 1 |         |   |  3   1          3   1   3
#  *     0---1---3---7     3---2  +-2-+          +-2-+-0-+-2-+
#  *         | 4 | 3 |                               3   1   3
#  *         0---2---6   Corners  Edges       Edges  +-2-+-0-+
#  */

# vtk configuration https://kitware.github.io/vtk-examples/site/Cxx/GeometricObjects/LinearCellDemo/
# /*
#  * vtk vertex and face orderning
#  *
#  *      7-------6       +   -   +       +-------+       +-------+
#  *     /|      /|      /|      /|      /   3   /        |   5   |
#  *    3-------2 |     + | -   + |     +-------+       +-------+ |
#  *    | |     | |     |0|     |1|                     | |     | |
#  *    | 4-----|-5     | +   - | +       +-------+     | +-----|-+
#  *    |/      |/      |/      |/       /   2   /      |   4   |
#  *    0-------1       +   -   +       +-------+       +-------+
#  *
#  *      Vertex                           Faces
#  */
def covert_mesh_poly_to_face(poly):
    # # hexaLab
    # faceConfig = [[0, 4, 6, 2],
    #               [1, 5, 7, 3],
    #               [0, 1, 5, 4],
    #               [2, 6, 7, 3],
    #               [0, 2, 3, 1],
    #               [4, 5, 7, 6]]
    # vtk
    faceConfig = [[0, 4, 5, 1],
                  [1, 5, 6, 2],
                  [0, 4, 7, 3],
                  [0, 1, 2, 3],
                  [6, 2, 3, 7],
                  [4, 5, 6, 7]]
    # if want to sync face vid order,
    # put small vid first,
    # and then only need to change direction of face to check sameness
    sortedFvids = []
    Faces = []
    allPolyFids = []
    for p in poly:
        polyFids = []
        for face in faceConfig:
            # get a face
            tmp = []
            for i in face:
                if p[i] in tmp:
                    continue
                tmp.append(p[i])
            # tmp = [p[i] for i in face]
            if len(tmp) < 3:
                # print(tmp)
                continue
            # try to get face id
            sortedTmp = sorted(tmp)
            try:
                fid = sortedFvids.index(sortedTmp)
                # polyFids.append(fid)
                if fid not in polyFids:
                    polyFids.append(fid)
            except:
                polyFids.append(len(Faces))
                Faces.append(tmp)
                sortedFvids.append(sortedTmp)
        allPolyFids.append(polyFids)
    # checking
    # for i, f in enumerate(Faces):
    #     for j, t in enumerate(Faces):
    #         if i == j:
    #             continue
    #         com = check_common_elements(f,t)
    #         if len(com) == len(f) or len(com) == len(t):
    #             print("--------------------- non conformal faces")
    #             print(com)
    #             print(f)
    #             print(t)

    return Faces, allPolyFids


def clean_vertices(vertices, faces):
    vid_use_map = [0] * len(vertices)
    vert_map = {}
    for face in faces:
        for v in face:
            vid_use_map[v] = 1
    for v in vid_use_map:
        if v == 0:
            print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")


# check whether the input face contains
# face 0,1,2,3 and also
# face 0,1,2
# faces = [[1]]
def check_sub_faces(faces):
    faces_stack = copy.deepcopy(faces)
    while faces_stack:
        face = faces_stack.pop()
        n = face[0]
        vids = face[1:]
        for otherF in faces_stack:
            on = otherF[0]
            ovids = otherF[1:]
            common_vids = check_common_elements(vids, ovids)
            lc = len(common_vids)
            if lc == n or lc == on:
                print("find a sub face pair ... ")
                print("face 0 : " + str(vids))
                print("face 1 : " + str(ovids))


# this function does not require face orientation and cell orientation
def create_basic_mesh(vertices, faces, cells, file_path="test.hybrid", mesh_type="basic_hexa_dominant"):
    new_cells = add_face_cell_orientation(cells)
    hm = BasicMesh(file_path, mesh_type)
    hm.build_mesh_from_vfc(vertices, faces, new_cells)
    return hm
    