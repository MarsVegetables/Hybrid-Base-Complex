from base_complex.base_complex import BaseComplex
from mesh_structure.multi_level_hybrid_mesh import MultiLevelHybridMesh


class BaseComplexToMeshConverter:
    def __init__(self, base_complex):
        self.base_complex = base_complex
        new_level = self.base_complex.mesh.mesh_level + 1
        new_file_path = self.base_complex.mesh.file_path
        new_mesh_type = self.base_complex.mesh.meshType
        self.new_mesh = MultiLevelHybridMesh(file_path=new_file_path,
                                             mesh_type=new_mesh_type,
                                             mesh_level=new_level)
        self.new_mesh.input_mesh = self.base_complex.mesh

    def covert(self):
        print("converting a base complex elements to mesh structure ...")
        # vertex
        for element in self.base_complex.componentVs:
            self.new_mesh.Vertices[element.element_id] = element
        # edge
        for element in self.base_complex.componentEs:
            self.new_mesh.Edges[element.element_id] = element
        # face
        for element in self.base_complex.componentFs:
            self.new_mesh.Faces[element.element_id] = element
        # cell
        for element in self.base_complex.componentCs:
            self.new_mesh.Cells[element.element_id] = element
        print("converting a base complex elements to mesh structure ... done")
        # build neighbor relationship
        self.new_mesh.build()
        # self.new_mesh.build_opposite_face_relationship()
        # self.new_mesh.build_parallel_face_relationship()
        # self.new_mesh.build_parallel_edge_relationship()

