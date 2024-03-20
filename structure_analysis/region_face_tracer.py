from mesh_structure.mesh_common import check_common_elements


class RegionFaceTracer:
    def __init__(self, mesh, avoid_eids, region_cids):
        self.mesh = mesh
        self.avoid_eids = avoid_eids
        self.region_cids = region_cids
        self.region_fids = set()
        self.inner_fids = set()
        self.face_patches = []
        # visited fids for tracing faces
        self.visited_fids = set()
        self.find_region_fids()
        self.find_inner_faces()
        self.get_face_patches()

    def is_face_valid(self, fid):
        face = self.mesh.Faces[fid]
        if check_common_elements(face.Eids, self.avoid_eids):
            return False
        return True

    def is_region_boundary_face(self, fid):
        face = self.mesh.Faces[fid]
        tmp = check_common_elements(face.neighboring_Cids,
                                    self.region_cids)
        # pass boundary fids
        if tmp == 1:
            return True
        return False

    def find_region_fids(self):
        # fids that not include any avoid eids
        for cid in self.region_cids:
            cell = self.mesh.Cells[cid]
            for fid in cell.Fids:
                if self.is_face_valid(fid):
                    self.region_fids.add(fid)

    # face layer - for self-parallel edge only
    # max matching will avoid face overlapping
    def find_inner_faces(self):
        self.inner_fids = set()
        for fid in self.region_fids:
            if self.is_region_boundary_face(fid):
                continue
            self.inner_fids.add(fid)
        return self.inner_fids

    def get_face_patches(self):
        self.visited_fids = set()
        self.face_patches = []
        for fid in self.region_fids:
            if fid in self.visited_fids:
                continue
            patch = self.trace_face_patch(fid)
            self.face_patches.append(patch)

    def trace_face_patch(self, fid):
        face_patch = set()
        face_stack = [fid]
        while face_stack:
            cur_fid = face_stack.pop()
            if cur_fid in self.visited_fids:
                continue
            self.visited_fids.add(cur_fid)
            face_patch.add(cur_fid)
            face = self.mesh.Faces[cur_fid]
            next_fids = check_common_elements(face.neighboring_Fids, self.region_fids)
            face_stack.extend(next_fids)
        return face_patch
