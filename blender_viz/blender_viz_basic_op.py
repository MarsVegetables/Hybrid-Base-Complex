# only for blender 3.0

import bpy
from mathutils import Matrix
# import bgl
import concurrent.futures
import math
from tqdm import tqdm
import colorsys
import time

def new_material(mid):
    mat = bpy.data.materials.new(name=mid)

    # mat = bpy.data.materials.get(id)

    # if mat is None:
    #     mat = bpy.data.materials.new(name=id)
    # mat.blend_method = "BLEND"
    mat.use_nodes = True
    mat.use_sss_translucency = True
    mat.use_screen_refraction = True
    mat.show_transparent_back = False
    # mat.use_backface_culling = True

    # if mat.node_tree:
    #     mat.node_tree.links.clear()
    #     mat.node_tree.nodes.clear()

    return mat


# https://www.katsbits.com/codex/alpha/
def new_shader(shader_id, type, rgba=[1., 0., 0., 1], diffusion=1, alpha=1):
    if (shader_id == "wireframe_color") or (shader_id == "camera_ui_color"):
        name = shader_id
    else:
        name = f"r_{rgba[0]:.3f}_g_{rgba[1]:.3f}_b_{rgba[2]:.3f}_a_{rgba[3]:.3f}"
    if name in bpy.data.materials.keys():
        mat = bpy.data.materials[name]
    else:
        mat = new_material(name)
    # mat.diffuse_color = rgba
    principled = mat.node_tree.nodes['Principled BSDF']
    principled.inputs["Subsurface Color"].default_value = rgba
    principled.inputs["Subsurface"].default_value = diffusion
    principled.inputs["Alpha"].default_value = alpha
    if alpha < 1:
        mat.blend_method = "HASHED"
        mat.shadow_method = "HASHED"
    return mat


def create_wireframe_geo_node(curve_radius=1):
    if "wireframe" in bpy.data.node_groups.keys():
        c_prof = bpy.data.node_groups["wireframe"].nodes["Curve Circle"]
        c_prof.inputs['Radius'].default_value = curve_radius
        return
    geo_group = bpy.data.node_groups.new("wireframe", "GeometryNodeTree")
    geo_group.inputs.new('NodeSocketGeometry', 'Geometry')
    geo_group.outputs.new('NodeSocketGeometry', 'Geometry')

    in_node = geo_group.nodes.new('NodeGroupInput')
    out_node = geo_group.nodes.new('NodeGroupOutput')

    mesh_to_c = geo_group.nodes.new("GeometryNodeMeshToCurve")
    c_to_mesh = geo_group.nodes.new("GeometryNodeCurveToMesh")
    c_prof = geo_group.nodes.new("GeometryNodeCurvePrimitiveCircle")

    c_prof.inputs['Radius'].default_value = curve_radius

    # wireframe color
    if "wireframe_color" not in bpy.data.materials.keys():
        mat = new_shader("wireframe_color", type="", rgba=[0.01, 0.01, 0.01, 1.])
    else:
        mat = bpy.data.materials["wireframe_color"]
    c_mat = geo_group.nodes.new("GeometryNodeSetMaterial")
    c_mat.inputs["Material"].default_value = mat

    geo_group.links.new(in_node.outputs['Geometry'], mesh_to_c.inputs["Mesh"])
    geo_group.links.new(mesh_to_c.outputs['Curve'], c_to_mesh.inputs["Curve"])
    geo_group.links.new(c_prof.outputs['Curve'], c_to_mesh.inputs["Profile Curve"])
    geo_group.links.new(c_to_mesh.outputs['Mesh'], c_mat.inputs['Geometry'])
    geo_group.links.new(c_mat.outputs['Geometry'], out_node.inputs['Geometry'])

    return geo_group.name


def make_hybrid_mesh_wireframe(vertices, edges, faces, collection_name='hybrid mesh'):
    # edges_t = []
    # faces_t = []

    # new mesh
    new_mesh = bpy.data.meshes.new('hybrid')
    new_mesh.from_pydata(vertices, edges, faces)
    # new_mesh.from_pydata(vertices, edges_t, faces)
    new_mesh.update()

    # make object from mesh
    obj_name = 'hybrid'
    new_object = bpy.data.objects.new(obj_name, new_mesh)

    geo_node = new_object.modifiers.new("show_wireframe", "NODES")
    geo_node.node_group = bpy.data.node_groups["wireframe"]

    # make collection
    new_collection = bpy.data.collections.new(collection_name)
    bpy.context.scene.collection.children.link(new_collection)

    # add object to scene collection
    new_collection.objects.link(new_object)
    # setObjGeometryToOrigin(obj_name)
    return obj_name


# add white background in blender
# https://www.youtube.com/watch?v=aegiN7XeLow
def add_materials_to_obj(ob, mat_name="color", rgba=(1, 0, 0, 1), alpha=2):
    if alpha > 1:
        alpha = rgba[-1]
    mat = new_shader(mat_name, type="diffuse", rgba=rgba, alpha=alpha)
    ob.active_material = mat


def make_poly_line(v_list, obj_name, collection_name=None, radius=0.1, rgba=None, circle=False):
    curve_data = bpy.data.curves.new(name=obj_name, type='CURVE')
    curve_data.dimensions = '3D'
    curve_data.bevel_depth = radius
    curve_data.bevel_mode = 'PROFILE'
    curve_data.use_fill_caps = True

    poly_line = curve_data.splines.new('POLY')
    poly_line.points.add(len(v_list) - 1)
    if circle:
        poly_line.use_cyclic_u = True
    for i, xyz in enumerate(v_list):
        poly_line.points[i].co = [xyz[0], xyz[1], xyz[2], 1]
        # poly_line.points[i].radius = radius

    obj_data = bpy.data.objects.new(obj_name, curve_data)
    # object origin
    # obj_data.bevel_depth = radius
    if rgba:
        add_materials_to_obj(obj_data, "line_color", rgba)

    if collection_name is not None:
        # make collection
        collection = get_blender_collection(collection_name)
        # bpy.context.scene.collection.objects.unlink(obj_data)
        collection.objects.link(obj_data)
    return poly_line


def create_cylinder_template():
    # Create a cylinder to use as a "template"
    bpy.ops.mesh.primitive_cylinder_add(
        radius=1,
        depth=1,
        enter_editmode=False,
    )
    cyl = bpy.context.object.data
    cyl.name = 'edge'
    bpy.context.object.hide_render = True
    bpy.data.objects.remove(bpy.context.object, do_unlink=True)
    return cyl


# add cylinder as an edge
def cylinder_between(v0, v1, r, collection=None, rgba=(1, 0, 0, 1)):
    x1 = v0[0]
    y1 = v0[1]
    z1 = v0[2]
    x2 = v1[0]
    y2 = v1[1]
    z2 = v1[2]
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    dist = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    bpy.ops.mesh.primitive_cylinder_add(
        radius=r,
        depth=dist,
        location=(dx / 2 + x1, dy / 2 + y1, dz / 2 + z1)
    )
    current_object = bpy.context.active_object

    add_materials_to_obj(current_object, "cylinder_color", rgba)

    # https://blender.stackexchange.com/questions/132112/whats-the-blender-2-8-command-for-adding-an-object-to-a-collection-using-python
    if collection is not None:
        bpy.context.scene.collection.objects.unlink(current_object)
        collection.objects.link(current_object)

    phi = math.atan2(dy, dx)
    theta = math.acos(dz / dist)

    current_object.rotation_euler[1] = theta
    current_object.rotation_euler[2] = phi

    return current_object


# https://blender.stackexchange.com/questions/233755/blender-taking-too-long-to-create-an-object-with-python-code
# create cylinder by using cylinder template
def cylinder_useTmp(v0, v1, r,
                    cylinder_tmp=None, collection=None,
                    rgba=(1, 0, 0, 1), edge_name = "edge"):
    if not cylinder_tmp:
        cylinder_between(v0, v1, r, collection=collection, rgba=rgba)
        return
    x1 = v0[0]
    y1 = v0[1]
    z1 = v0[2]
    x2 = v1[0]
    y2 = v1[1]
    z2 = v1[2]
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    dist = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    # cylinder_tmp = bpy.ops.mesh.primitive_cylinder_add(
    #     radius = r,
    #     depth = dist,
    #     location = (dx/2 + x1, dy/2 + y1, dz/2 + z1)
    # )

    cyl_mesh = cylinder_tmp.copy()
    cyl_mesh.transform(Matrix.Diagonal([r, r, dist, 1]))
    cyl_mesh.name = edge_name

    cyl_obj = bpy.data.objects.new(cyl_mesh.name, cyl_mesh)
    cyl_obj.hide_render = False
    add_materials_to_obj(cyl_obj, "cylinder_color", rgba)

    # https://blender.stackexchange.com/questions/132112/whats-the-blender-2-8-command-for-adding-an-object-to-a-collection-using-python
    if collection is not None:
        collection.objects.link(cyl_obj)

    phi = math.atan2(dy, dx)
    if dist > 0:
        theta = math.acos(dz / dist)
    else:
        theta = 0

    cyl_obj.location = (dx / 2 + x1, dy / 2 + y1, dz / 2 + z1)
    cyl_obj.rotation_euler[1] = theta
    cyl_obj.rotation_euler[2] = phi

    return cyl_obj


def get_blender_collection(collection_name='collection'):
    new_collection = bpy.data.collections.get(collection_name)
    if new_collection is None:
        # make collection
        new_collection = bpy.data.collections.new(collection_name)
    # if it is not linked to scene colleciton treelink it
    if not bpy.context.scene.user_of_id(new_collection):
        bpy.context.scene.collection.children.link(new_collection)
    return new_collection


def make_edges(vertices, edges, r=0.02,
               name='non hexa mesh', rgba=(1,0,0,1), edge_names=[],
               edge_template=None):
    # print("making edges in blender ...")
    # make collection
    new_collection = get_blender_collection(name)
    if edge_template:
        cyl_tmp = edge_template
    else:
        cyl_tmp = create_cylinder_template()
    # add object to scene collection
    use_name = False
    if edge_names:
        use_name = True
    for i, e in enumerate(edges):
        if e[0] == e[1]:
            continue
        v0 = vertices[e[0]]
        v1 = vertices[e[1]]
        # cylinder_between(v0, v1, r, collection = new_collection, rgba = rgba)
        if use_name:
            edge_name = "edge_{}".format(edge_names[i])
            cyl_obj = cylinder_useTmp(v0, v1, r=r,
                                      cylinder_tmp=cyl_tmp,
                                      collection=new_collection,
                                      rgba=rgba,
                                      edge_name=edge_name)
        else:
            cyl_obj = cylinder_useTmp(v0, v1, r=r, cylinder_tmp=cyl_tmp,
                                      collection=new_collection, rgba=rgba)
    # bpy.context.scene.collection.objects.unlink(current_object)
    # collection.objects.link(current_object)


# face_vid_groups contain many faces' vids
# n * n array
# createFacesInOneMeshWithSameColor
def visualize_faces(vertices, face_vid_groups, rgba=None,
                    face_name='face', collection_name='faces'):
    edges_t = []

    # make collection
    new_collection = get_blender_collection(collection_name)

    # new mesh
    new_mesh = bpy.data.meshes.new(face_name)
    new_mesh.from_pydata(vertices, edges_t, face_vid_groups)
    # new_mesh.from_pydata(vertices, edges_t, faces)
    new_mesh.update()

    # make object from mesh
    new_object = bpy.data.objects.new(face_name, new_mesh)
    if rgba:
        add_materials_to_obj(new_object, "group_color", rgba)

    # add object to scene collection
    new_collection.objects.link(new_object)


# hue is 0 to 1
def color_id_to_rgba(color_id=0, total_color=8, alpha=1):
    hue = color_id / total_color
    if hue > 1:
        hue = 1
    rgb = colorsys.hsv_to_rgb(hue, 1, 1)
    return [rgb[0], rgb[1], rgb[2], alpha]


def create_sphere_tmp(radius=0.1):
    bpy.ops.mesh.primitive_ico_sphere_add(radius=radius, location=[0,0,0])
    current_object = bpy.context.active_object
    current_object.hide_render = True
    return current_object.data


def make_sphere(xyz, radius=0.1, sphere_name='vert',
                collection_name='spheres', rgba=(1, 0, 0, 1),
                sphere_template=None, size_weight=1):
    # make collection
    new_collection = get_blender_collection(collection_name)

    if not sphere_template:
        sph_mesh = create_sphere_tmp(radius)
    else:
        sph_mesh = sphere_template.copy()
    sph_obj = bpy.data.objects.new(sphere_name, sph_mesh)
    # sph_obj.hide_render = False
    # bpy.data.objects.remove(sph_mesh, do_unlink=True)
    sph_obj.location = xyz
    sph_obj.scale = size_weight * sph_obj.scale
    add_materials_to_obj(sph_obj, mat_name="sphere_color", rgba=rgba)
    # bpy.context.scene.collection.objects.unlink(sph_obj)
    # add object to scene collection
    new_collection.objects.link(sph_obj)


# set world data
def set_world_color(color=(1, 1, 1, 1)):
    bpy.data.worlds['World'].node_tree.nodes['Background'].inputs['Color'].default_value = color
    bpy.data.worlds['World'].node_tree.nodes['Background'].inputs['Strength'].default_value = 0.1


def set_color_management(view_transform='Standard', render_square_resolution=True):
    scn = bpy.context.scene
    # ('Standard', 'Filmic', 'Filmic Log', 'Raw', 'False Color')
    scn.view_settings.view_transform = view_transform
    scn.render.film_transparent = True
    scn.render.image_settings.file_format = 'PNG'
    scn.render.image_settings.color_mode = 'RGBA'
    scn.render.image_settings.color_depth = '16'
    scn.render.image_settings.compression = 15
    scn.display_settings.display_device = 'sRGB'
    scn.sequencer_colorspace_settings.name = 'sRGB'

    scn.render.engine = 'BLENDER_EEVEE'
    scn.eevee.use_gtao = True
    # scn.eevee.taa_render_samples = 1280
    scn.eevee.taa_render_samples = 128
    if render_square_resolution:
        scn.render.resolution_x = 1080
        scn.render.resolution_y = 1080
    else:
        scn.render.resolution_x = 1920
        scn.render.resolution_y = 1080
    # bpy.context.scene.render.cycles.device = 'GPU'

    # enable ao
    scn.world.light_settings.use_ambient_occlusion = True  # turn AO on
    scn.world.light_settings.ao_factor = 0.5  # set it to 0.5


def remove_all_materials():
    default_materials = ["Dots Stroke", "Material", "wireframe_color", "camera_ui_color"]
    all_materials = []
    for material in bpy.data.materials:
        if material.name_full in default_materials:
            continue
        # https://blender.stackexchange.com/questions/260190/id-user-decrement-error-raised-when-i-try-to-clear-users-from-a-screen-what-doe
        # material.user_clear()
        # bpy.data.materials.remove(material)
        all_materials.append(material)
    bpy.data.batch_remove(all_materials)


def remove_blender_mesh(blender_meshes):
    for obj in blender_meshes:
        # Delete the meshes
        bpy.data.meshes.remove(obj.data)
def batch_mesh_remove(blender_meshes):
    n = len(blender_meshes)
    ideal_batch_size = 300
    worker_num = math.ceil(n/ideal_batch_size)
    batch_size = int(n / worker_num)
    if worker_num <= 5:
        remove_blender_mesh(blender_meshes)
        return
    print(f"batch removing blender meshes ...")
    print(f"number of worker : {worker_num}, batch_size : {batch_size}")
    pool = concurrent.futures.ThreadPoolExecutor(max_workers=worker_num)
    for i in range(worker_num):
        start_i = i * batch_size
        if i == worker_num - 1:
            end_i = len(blender_meshes)
        else:
            end_i = (i + 1) * batch_size
        batch_meshes = blender_meshes[start_i:end_i]
        pool.submit(remove_blender_mesh, batch_meshes)
    pool.shutdown(wait=True)
    print("batch removing blender meshes ... done")

def remove_collection(collection_name, remove_meshes = False):
    time_s = time.time()
    if not collection_name:
        return
    if collection_name not in bpy.data.collections:
        return
    # Get the collection from its name
    collection = bpy.data.collections[collection_name]
    if not collection:
        return

    if remove_meshes:
        meshes = []
        obj_ids = []
        for obj in collection.objects:
            # Delete the object
            if obj.type == 'MESH':
                meshes.append(obj)
            else:
                obj_ids.append(obj)
        # bpy.data.meshes.batch_remove(mesh_ids)
        # for obj in meshes:
        #     # Delete the meshes
        #     bpy.data.meshes.remove(obj.data)
        batch_mesh_remove(meshes)
        bpy.data.batch_remove(obj_ids)
    else:
        scene = bpy.data.scenes['Scene']
        scene.collection.children.unlink(collection)
    bpy.data.collections.remove(collection)
    time_e = time.time()
    print("collection remove take {:.2f} second".format(time_e - time_s))


def collection_render_control(collection_name, hide_render=False):
    if collection_name not in bpy.data.collections:
        return hide_render
    coll = bpy.data.collections[collection_name]
    coll.hide_render = hide_render
    # for c in coll.objects:
    #     c.hide_render = hide_render
    return hide_render
