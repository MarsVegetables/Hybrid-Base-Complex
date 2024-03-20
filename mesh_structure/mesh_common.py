'''
Author: Lei Si
Date: 2022-08-22 11:48:04
LastEditTime: 2022-09-09 01:29:04
LastEditors: Lei Si
Description: 
FilePath: \hexa_dominant_layer_simplification\meshLoader\meshCommon.py
YooooHooooo
'''
import numpy as np
from tqdm import tqdm
# gc : https://stackoverflow.com/questions/2473783/is-there-a-way-to-circumvent-python-list-append-becoming-progressively-slower
import gc
from copy import deepcopy


# faces to edges 
def face_list_to_edge_list(faces_vids):
    # vertex number is edge number in a face list
    print("converting faces to edges ...")
    all_edges = []
    faces_eids = []
    for f in tqdm(faces_vids):
        gc.disable()
        face_eids = []
        f_edges = face_to_edges(f)
        for edge in f_edges:
            try:
                eid = all_edges.index(edge)
            except ValueError:
                eid = len(all_edges)
                all_edges.append(edge)
            face_eids.append(eid)
        faces_eids.append(face_eids)
        gc.enable()
    print("converting faces to edges ... done")
    return all_edges, faces_eids


# split a n-vertices face to n edges
def face_to_edges(face_vids):
    edges = []
    n = len(face_vids)
    vids = face_vids
    for i in range(n):
        v0 = i % n
        v1 = (i + 1) % n
        edge = [vids[v0], vids[v1]]
        edge.sort()
        edges.append(edge)
    return edges


def check_common_elements(list_0, list_2):
    return set(list_0) & set(list_2)


# mixed prodcut can be considered as the volume of a hexa
# all vertex point to v0
def mixed_product(v0, v1, v2, v3):
    a = normalize(np.array(v1) - np.array(v0))
    b = normalize(np.array(v2) - np.array(v0))
    c = normalize(np.array(v3) - np.array(v0))
    bc_cross = np.cross(b, c)
    abc_mixed = np.dot(a, bc_cross)
    return abc_mixed


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


def add_face_cell_orientation(cells):
    # each cell contain face id for the cell
    # add Face Orientation And Cell orientation
    new_cells = []
    for cell in cells:
        new_c = list()
        new_c.append(cell)
        new_c.append([0 for i in range(len(cell))]) # face orientation
        new_c.append(0) # cell orientation
        new_cells.append(new_c)
    return new_cells


def reduce_id_groups(ll):
    l = deepcopy(ll)
    s = set()
    for i in l:
        s.update(i)
    for i in s:
        components = [x for x in l if i in x]
        for j in components:
            l.remove(j)
        new_group = set()
        for c in components:
            new_group.update(c)
        l += [list(new_group)]
    return l


def varify_and_clean_file_data(faces, polys):
    new_faces = []
    for f in faces:
        n = f[0]
        vids = f[1:]
        if n != len(vids):
            print("not match face ----- ")
            print(f)
        new_faces.append(vids)

    new_polys = []
    for p in polys:
        new_p = [p[0][1:], p[1][1:], p[2]]
        for i in range(len(p) - 1):
            if p[i][0] != len(new_p[i]):
                print("not match poly ----")
                print(p)
        new_polys.append(new_p)
    return new_faces, new_polys


def remove_duplicate_faces_and_remapping(faces, polys):
    v_to_f = {}
    for i, face_vids in enumerate(faces):
        for vid in face_vids:
            if vid not in v_to_f:
                v_to_f[vid] = set()
            v_to_f[vid].add(i)
    face_new_id = [-1] * len(faces)
    new_face_list = []
    for i, face_vids in enumerate(faces):
        # marked faces
        if face_new_id[i] != -1:
            continue
        new_id = len(new_face_list)
        face_new_id[i] = new_id
        new_face_list.append(face_vids)
        # search faces
        vid = face_vids[0]
        fids = v_to_f[vid]
        for fid in fids:
            if set(face_vids) == set(faces[fid]):
                face_new_id[fid] = new_id
    # remap cells
    new_polys = []
    for poly in polys:
        new_fids = set()
        for fid in poly[0]:
            new_fids.add(face_new_id[fid])
        if len(new_fids) > 3:
            poly[0] = list(new_fids)
            poly[1] = [0] * len(new_fids)
            new_polys.append(poly)
    return new_face_list, new_polys


def find_duplicate_elements(element_list, order_matter=True):
    if order_matter:
        new_list = element_list
    else:
        new_list = [set(x) for x in element_list]
    return len(element_list) - len(np.unique(new_list, axis=0))
    # existing_elements = []
    # duplicate_count = 0
    # for e in element_list:
    #     ee = set(e)
    #     if order_matter:
    #         ee = e
    #     if ee in existing_elements:
    #         duplicate_count += 1
    #     else:
    #         existing_elements.append(ee)
    # return duplicate_count









