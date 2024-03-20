import bpy

import math
import mathutils
from pathlib import Path
import os

from blender_viz.blender_viz_basic_op import (get_blender_collection,
                                              set_world_color, set_color_management,
                                              remove_collection, remove_all_materials,
                                              new_shader)

# https://blenderartists.org/t/how-to-create-a-mathutils-matrix-from-location-rotation-scale/1242751/2
# https://blender.stackexchange.com/questions/254903/align-active-camera-to-view-using-python


class CameraCapture:
    def __init__(self, file_path, file_name, show_mesh_info=True):
        self.file_path = file_path
        self.file_name = file_name
        self.camera = None
        self.camera_collection = None
        self.obj_to_camera = None
        self.target_obj = None
        self.light = None

        # camera UI
        self.show_info = show_mesh_info
        self.camera_ui_node_group = None
        self.string_obj = None
        self.string_node = None
        #
        self.image_idx = 0
        set_world_color()
        set_color_management()

    def update_camera_and_light(self, target_name="Cube"):
        self.create_camera()
        self.add_light()
        self.move_view_to_object(target_name)
        if self.show_info:
            self.create_camera_ui()

    def create_camera(self, name="Camera", collection_name="World"):
        self.camera = bpy.context.scene.objects.get(name)
        if self.camera:
            self.camera_collection = get_blender_collection(collection_name)
            return self.camera
        self.camera = bpy.data.objects.new(name, bpy.data.cameras.new(name))

        # make collection
        new_collection = get_blender_collection(collection_name)

        # add object to scene collection
        new_collection.objects.link(self.camera)
        self.camera_collection = new_collection
        return self.camera

    def move_view_to_object(self, target_name="hybrid"):
        bpy.ops.object.select_all(action='DESELECT')
        self.target_obj = bpy.context.scene.objects[target_name]
        self.target_obj.select_set(True)
        bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='BOUNDS')
        # set camera
        bpy.context.scene.camera = self.camera
        # must update rotation first
        self.adjust_camera_rotation()
        self.adjust_camera_position()
        self.update_light_direction()

    def adjust_camera_position(self):
        if not self.target_obj:
            return
        # change view to the object
        bpy.ops.view3d.camera_to_view_selected()
        obj_to_camera = self.camera.location - self.target_obj.location
        self.obj_to_camera = obj_to_camera
        new_camera_position = (1.1 * obj_to_camera) + self.target_obj.location
        self.camera.location = new_camera_position

    def adjust_camera_rotation(self):
        if not self.target_obj:
            return
        rot_quat = self.target_obj.matrix_world.to_4x4()
        # rotx = mathutils.Matrix.Rotation(math.radians(-30), 4, 'X')
        # roty = mathutils.Matrix.Rotation(math.radians(20), 4, 'Y')
        # rotz = mathutils.Matrix.Rotation(math.radians(0), 4, 'Z')
        # camera_rot = rotz @ roty @ rotx @ rot_quat
        camera_rotx = mathutils.Matrix.Rotation(math.radians(60), 4, 'X')
        camera_rotz = mathutils.Matrix.Rotation(math.radians(45), 4, 'Z')
        camera_rot = camera_rotz @ camera_rotx @ rot_quat
        self.camera.rotation_euler = camera_rot.to_euler()

    def add_light(self, name="Light", collection_name="World"):
        self.light = bpy.context.scene.objects.get(name)
        if self.light:
            if self.light.data.energy != 5.0:
                self.light.data.energy = 5.0
        else:
            light_data = bpy.data.lights.new(name=name, type='SUN')
            light_data.energy = 5.0
            light_data.use_shadow = False
            new_collection = get_blender_collection(collection_name)
            self.light = bpy.data.objects.new(name=name, object_data=light_data)
            new_collection.objects.link(self.light)
        return self.light

    def update_light_direction(self):
        if self.light and self.target_obj:
            rot_quat = self.target_obj.matrix_world.to_4x4()
            # rotx = mathutils.Matrix.Rotation(math.radians(-30), 4, 'X')
            # roty = mathutils.Matrix.Rotation(math.radians(20), 4, 'Y')
            # rotz = mathutils.Matrix.Rotation(math.radians(0), 4, 'Z')
            # rot = rotz @ roty @ rotx @ rot_quat
            rotx = mathutils.Matrix.Rotation(math.radians(60), 4, 'X')
            rotz = mathutils.Matrix.Rotation(math.radians(45), 4, 'Z')
            rot = rotz @ rotx @ rot_quat
            self.light.rotation_euler = rot.to_euler()

    def render_view(self, img_name='', target_info=None, index_name=False):
        file_folder_name = self.file_name.split(".")[0]
        # file_folder = Path(file_path, f'/renders/{img_name}')
        file_folder = Path(self.file_path)
        file_folder = file_folder.joinpath('renders')
        if not os.path.exists(file_folder):
            os.mkdir(file_folder)
        file_folder = file_folder.joinpath(file_folder_name)
        if not os.path.exists(file_folder):
            os.mkdir(file_folder)
        if index_name or (not img_name):
            bpy.context.scene.render.filepath = str(file_folder.joinpath(str(self.image_idx)))
            self.image_idx += 1
        else:
            bpy.context.scene.render.filepath = str(file_folder.joinpath(img_name))
        # update mesh information
        if target_info and self.show_info:
            all_info = target_info
            if img_name:
                all_info += img_name
            self.update_string(all_info)
        # Render Scene and store the scene
        bpy.ops.render.render(write_still=True)

    def clean_old_collection(self, reset_img_idx=False):
        if reset_img_idx:
            self.image_idx = 0
        for co in bpy.data.collections:
            if co.name_full == self.camera_collection.name_full:
                continue
            remove_collection(co.name_full, remove_meshes=True)
        remove_all_materials()

    def create_camera_ui(self, node_group_name="camera_ui"):
        # node_group
        if node_group_name in bpy.data.node_groups.keys():
            return
        geo_group = bpy.data.node_groups.new(node_group_name, "GeometryNodeTree")
        geo_group.inputs.new('NodeSocketGeometry', 'Geometry')
        geo_group.outputs.new('NodeSocketGeometry', 'Geometry')
        nodes = geo_group.nodes

        in_node = nodes.new('NodeGroupInput')
        out_node = nodes.new('NodeGroupOutput')

        string_node = nodes.new('FunctionNodeInputString')
        self.string_node = string_node
        string_2_curve = nodes.new('GeometryNodeStringToCurves')
        transform_geo = nodes.new('GeometryNodeTransform')
        fill_curve = nodes.new('GeometryNodeFillCurve')

        # text color
        mat = new_shader("camera_ui_color", type="", rgba=(0.0, 0.0, 0.0, 1))
        c_mat = geo_group.nodes.new("GeometryNodeSetMaterial")
        c_mat.inputs["Material"].default_value = mat

        geo_group.links.new(string_node.outputs['String'], string_2_curve.inputs["String"])
        geo_group.links.new(string_2_curve.outputs['Curve Instances'], transform_geo.inputs["Geometry"])
        geo_group.links.new(transform_geo.outputs['Geometry'], fill_curve.inputs['Curve'])
        geo_group.links.new(fill_curve.outputs['Mesh'], c_mat.inputs['Geometry'])
        geo_group.links.new(c_mat.outputs['Geometry'], out_node.inputs['Geometry'])

        self.camera_ui_node_group = geo_group

        # add to scene
        new_mesh = bpy.data.meshes.new(node_group_name)
        new_obj = bpy.data.objects.new(node_group_name, new_mesh)
        self.string_obj = new_obj

        # add text to obj
        geo_node = new_obj.modifiers.new(node_group_name, "NODES")
        geo_node.node_group = bpy.data.node_groups[node_group_name]

        self.camera_collection.objects.link(new_obj)
        # rotate text
        new_obj.rotation_euler = self.camera.rotation_euler

        # set text position
        cam_angle_x = self.camera.data.angle_x
        cam_angle_y = self.camera.data.angle_y
        # local t_x and t_y
        t_z = 0.5 * self.obj_to_camera.length
        t_x = t_z * math.tan(0.5 * cam_angle_x)
        t_y = t_z * math.tan(0.5 * cam_angle_y)
        local_xyz = mathutils.Vector((-0.98 * t_x, 0.75 * t_y, -t_z))
        cam_matrix = self.camera.matrix_world.to_quaternion().to_matrix().to_4x4()
        world_xyz = cam_matrix @ local_xyz

        new_obj.location = self.camera.location + world_xyz

        # scale
        s_x = (0.06 * t_x) / new_obj.dimensions.x
        new_obj.scale = (s_x, s_x, s_x)

    def update_string(self, mesh_info):
        self.string_node.string = mesh_info


