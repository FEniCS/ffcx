"This module provides functionality for plotting finite elements."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2010-12-07"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

__all__ = ["plot"]

from numpy import dot, cross, array, sin, cos, pi
from numpy.linalg import norm

from ffc.fiatinterface import create_element
from ffc.log import warning, error

# Import Soya3D
try:
    import soya
    from soya.sphere import Sphere
    _soya_imported = True
except:
    _soya_imported = False

def plot(element):
    "Plot finite element."

    # Check if Soya3D has been imported
    if not _soya_imported:
        warning("Unable to plot element, Soya3D not available (install package python-soya).")
        return

    # Create cell model
    cell = create_cell_model(element)

    # Create dof models
    dofs = create_dof_models(element)

    # Render plot window
    render([cell] + dofs)

def render(models):
    "Render given list of models."

    # Note that we view from the positive z-axis, and not from the
    # negative y-axis. This should make no difference since the
    # element dofs are symmetric anyway and it plays better with
    # the default camera settings in Soya.

    # Initialize Soya
    soya.init("FFC plot", sound=1)

    # Create scene
    scene = soya.World()
    scene.atmosphere = soya.Atmosphere()
    scene.atmosphere.bg_color = (1.0, 1.0, 1.0, 1.0)

    # Define rotation around y-axis
    class RotatingBody(soya.Body):
        def advance_time(self, proportion):
            self.rotate_y(2.0 * proportion)

    # Add all models
    for model in models:
        rotating = RotatingBody(scene, model)

    # Set light
    light = soya.Light(scene)
    light.set_xyz(2.0, 5.0, -1.0)
    light.cast_shadow = 1
    light.shadow_color = (0.0, 0.0, 0.0, 0.5)

    # Set camera
    camera = soya.Camera(scene)
    camera.set_xyz(0.0, 10, 50.0)
    camera.ortho = 0
    camera.fov = 1.8
    p = camera.position()
    p.set_xyz(0.0, 0.4, 0.0)
    camera.look_at(p)
    soya.set_root_widget(camera)

    # Main loop
    soya.MainLoop(scene).main_loop()

def tangents(n):
    "Return normalized tangent vectors for plane defined by given vector."

    # Figure out which vector to take cross product with
    eps = 1e-10
    e = array((1.0, 0.0, 0.0))
    if norm(cross(n, e)) < eps:
        e = array((0.0, 1.0, 0.0))

    # Take cross products and normalize
    t0 = cross(n, e)
    t0 = t0 / norm(t0)
    t1 = cross(n, t0)
    t1 = t1 / norm(t0)

    return t0, t1

def Cylinder(scene, p0, p1, r, color=(0.0, 0.0, 0.0, 1.0)):
    "Return model for cylinder from p0 to p1 with radius r."

    # Convert to NumPy array
    if isinstance(p0, soya.Vertex):
        p0 = array((p0.x, p0.y, p0.z))
        p1 = array((p1.x, p1.y, p1.z))
    else:
        p0 = array(p0)
        p1 = array(p1)

    # Get tangent vectors for plane
    n = p0 - p1
    n = n / norm(n)
    t0, t1 = tangents(n)

    # Traverse the circles
    num_steps = 10
    dtheta = 2.0*pi / float(num_steps)
    for i in range(num_steps):

        # Compute coordinates for square
        dx0 = cos(i*dtheta)*t0 + sin(i*dtheta)*t1
        dx1 = cos((i + 1)*dtheta)*t0 + sin((i + 1)*dtheta)*t1
        x0 = p0 + r*dx0
        x1 = p0 + r*dx1
        x2 = p1 + r*dx0
        x3 = p1 + r*dx1

        # Cover square by two triangles
        v0 = soya.Vertex(scene, x0[0], x0[1], x0[2], diffuse=color)
        v1 = soya.Vertex(scene, x1[0], x1[1], x1[2], diffuse=color)
        v2 = soya.Vertex(scene, x2[0], x2[1], x2[2], diffuse=color)
        v3 = soya.Vertex(scene, x3[0], x3[1], x3[2], diffuse=color)
        f0 = soya.Face(scene, (v0, v1, v2))
        f1 = soya.Face(scene, (v1, v2, v3))
        f0.double_sided = 1
        f1.double_sided = 1

    # Extract model
    model = scene.to_model()

    return model

def Arrow(scene, x, n, l=0.3):
    "Return model for arrow from x in direction n."

    # Convert to Numpy arrays
    x = array(x)
    n = array(n)

    # Create cylinders
    l1 = Cylinder(scene, x, x + l*n, 0.04*l)

    # Extract model
    return scene.to_model()

def UnitTetrahedron(color=(0.0, 1.0, 0.0, 0.5)):
    "Return model for unit tetrahedron."

    # Create separate scene (since we will extract a model, not render)
    scene = soya.World()

    # Create vertices
    v0 = soya.Vertex(scene, 0.0, 0.0, 0.0, diffuse=color)
    v1 = soya.Vertex(scene, 1.0, 0.0, 0.0, diffuse=color)
    v2 = soya.Vertex(scene, 0.0, 1.0, 0.0, diffuse=color)
    v3 = soya.Vertex(scene, 0.0, 0.0, 1.0, diffuse=color)

    # Create edges
    e0 = Cylinder(scene, v0, v1, 0.005)
    e0 = Cylinder(scene, v0, v2, 0.005)
    e0 = Cylinder(scene, v0, v3, 0.005)
    e0 = Cylinder(scene, v1, v2, 0.005)
    e0 = Cylinder(scene, v1, v3, 0.005)
    e0 = Cylinder(scene, v2, v3, 0.005)

    # Create faces
    f0 = soya.Face(scene, (v1, v2, v3))
    f1 = soya.Face(scene, (v0, v2, v3))
    f2 = soya.Face(scene, (v0, v1, v3))
    f3 = soya.Face(scene, (v0, v1, v2))

    # Make faces double sided
    f0.double_sided = 1
    f1.double_sided = 1
    f2.double_sided = 1
    f3.double_sided = 1

    # Extract model
    model = scene.to_model()

    return model

def UnitTriangle(color=(0.0, 1.0, 0.0, 0.5)):
    "Return model for unit tetrahedron."

    # Create separate scene (since we will extract a model, not render)
    scene = soya.World()

    # Create vertice
    v0 = soya.Vertex(scene, 0.0, 0.0, 0.0, diffuse=color)
    v1 = soya.Vertex(scene, 1.0, 0.0, 0.0, diffuse=color)
    v2 = soya.Vertex(scene, 0.0, 1.0, 0.0, diffuse=color)

    # Create face
    f = soya.Face(scene, (v0, v1, v2))

    # Make face double sided
    f.double_sided = 1

    # Extract model
    model = scene.to_model()

    return model

def NormalEvaluation(x, n):
    "Return model for normal evaluation at given point in given direction."

    # Create separate scene (since we will extract a model, not render)
    scene = soya.World()

    # Normalize
    n = array(n)
    n = 0.75 * n / norm(n)

    # Flip normal if necessary
    if not pointing_outwards(x, n):
        n = -n

    # Create arrow
    arrow = Arrow(scene, x, n)

    # Extract model
    model = scene.to_model()

    return model

def PointEvaluation(x):
    "Return model for point evaluation at given point."

    # Create separate scene (since we will extract a model, not render)
    scene = soya.World()

    # Define material (color) for the sphere
    material = soya.Material()
    material.diffuse = (0.0, 0.0, 0.0, 1.0)

    # Create sphere
    sphere = Sphere(scene, material=material)

    # Scale and moveand move to coordinate
    sphere.scale(0.05, 0.05, 0.05)
    p = sphere.position()
    p.set_xyz(x[0], x[1], x[2])
    sphere.move(p)

    # Extract model
    model = scene.to_model()

    return model

def create_cell_model(element):
    "Create Soya3D model for cell."

    # Create model based on domain type
    domain = element.cell().domain()
    if domain == "triangle":
        return UnitTriangle()
    elif domain == "tetrahedron":
        return UnitTetrahedron()

    error("Unable to plot element, unhandled cell type: %s" % str(domain))

def create_dof_models(element):
    "Create Soya3D models for dofs."

    # Create FIAT element and get dofs
    fiat_element = create_element(element)
    dofs = fiat_element.dual_basis()

    # Iterate over dofs and add models
    models = []
    for dof in dofs:

        # Get type of dof
        dof_type = dof.get_type_tag()

        # Create model based on dof type
        L = dof.get_point_dict()
        if dof_type == "PointEval":

            # Point evaluation, just get point
            points = L.keys()
            if not len(points) == 1:
                error("Strange dof, single point expected for point evaluation.")
            x = points[0]

            # Generate model
            models.append(PointEvaluation(x))

        elif dof_type == "PointScaledNormalEval":

            # Normal evaluation, get point and normal
            points = L.keys()
            if not len(points) == 1:
                error("Strange dof, single point expected for point evaluation.")
            x = points[0]
            n = [xx[0] for xx in L[x]]

            # Generate model
            models.append(NormalEvaluation(x, n))

        elif dof_type == "FrobeniusIntegralMoment":

            warning("Not plotting interior moment for now.")

        else:
            error("Unable to plot dof, unhandled dof type: %s" % str(dof_type))

    return models

def pointing_outwards(x, n):
    "Check if n is pointing inwards, used for flipping dofs."
    eps = 1e-10
    x = array(x) + 0.1*array(n)
    return x[0] < -eps or x[1] < -eps or x[2] < -eps or x[2] > 1.0 - x[0] - x[1] + eps

