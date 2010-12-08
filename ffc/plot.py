"This module provides functionality for plotting finite elements."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2010-12-07"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-12-08

__all__ = ["plot"]

from numpy import dot, cross, array, sin, cos, pi, sqrt
from numpy.linalg import norm

from ffc.fiatinterface import create_element
from ffc.log import warning, error, info

# Import Soya3D
try:
    import soya
    from soya.sphere import Sphere
    from soya.label3d import Label3D
    _soya_imported = True
except:
    _soya_imported = False

# Colors for elements
element_colors = {"Argyris":                  (0.45, 0.70, 0.80),
                  "Arnold-Winther":           (0.00, 0.00, 1.00),
                  "Brezzi-Douglas-Marini":    (1.00, 1.00, 0.00),
                  "Crouzeix-Raviart":         (1.00, 0.25, 0.25),
                  "Discontinuous Lagrange":   (0.00, 0.25, 0.00),
                  "Hermite":                  (0.50, 1.00, 0.50),
                  "Lagrange":                 (0.00, 1.00, 0.00),
                  "Mardal-Tai-Winther":       (1.00, 0.10, 0.90),
                  "Morley":                   (0.40, 0.40, 0.40),
                  "Nedelec 1st kind H(curl)": (0.90, 0.30, 0.00),
                  "Raviart-Thomas":           (0.90, 0.60, 0.00)}

def plot(element):
    "Plot finite element."

    # Check if Soya3D has been imported
    if not _soya_imported:
        warning("Unable to plot element, Soya3D not available (install package python-soya).")
        return

    # Create cell model
    cell, is3d = create_cell_model(element)

    # Create dof models
    dofs, num_moments = create_dof_models(element)

    # Render plot window
    render(element, [cell] + dofs, num_moments, is3d)

def render(element, models, num_moments, is3d=True):
    "Render given list of models."

    # Note that we view from the positive z-axis, and not from the
    # negative y-axis. This should make no difference since the
    # element dofs are symmetric anyway and it plays better with
    # the default camera settings in Soya.

    # Initialize Soya
    if element.degree() is not None:
        title = "%s degree %d on a %s" % (element.family(), element.degree(), element.cell().domain())
    else:
        title = "%s degree on a %s" % (element.family(), element.cell().domain())
    soya.init(title)

    # Create scene
    scene = soya.World()
    scene.atmosphere = soya.Atmosphere()
    scene.atmosphere.bg_color = (1.0, 1.0, 1.0, 1.0)

    # Not used, need to manually handle rotation
    #label = Label3D(scene, text=str(num_moments), size=0.005)
    #label.set_xyz(1.0, 1.0, 1.0)
    #label.set_color((0.0, 0.0, 0.0, 1.0))

    # Define rotation around y-axis
    if is3d:
        class RotatingBody(soya.Body):
            def advance_time(self, proportion):
                self.rotate_y(2.0 * proportion)
    else:
        class RotatingBody(soya.Body):
            def advance_time(self, proportion):
                self.rotate_z(2.0 * proportion)

    # Add all models
    for model in models:
        rotating = RotatingBody(scene, model)

    # Set light
    light = soya.Light(scene)
    if is3d:
        light.set_xyz(1.0, 5.0, 5.0)
    else:
        light.set_xyz(0.0, 0.0, 1.0)
    light.cast_shadow = 1
    light.shadow_color = (0.0, 0.0, 0.0, 0.5)

    # Set camera
    camera = soya.Camera(scene)
    camera.ortho = 0
    camera.set_xyz(0.0, 10, 50.0)
    p = camera.position()
    if is3d:
        camera.fov = 2.1
        p.set_xyz(0.0, 0.4, 0.0)
    else:
        camera.fov = 2.6
        p.set_xyz(0.0, 0.0, 0.0)
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

def Arrow(scene, x, n, center=False):
    "Return model for arrow from x in direction n."

    # Convert to Numpy arrays
    x = array(x)
    n = array(n)

    # Get tangents
    t0, t1 = tangents(n)

    # Dimensions for arrow
    l = 0.3
    r = 0.04*l
    R1 = 0.1*l
    R2 = 0.2*l

    # Create cylinders
    x0 = x
    x1 = x + l*n
    x2 = x1 - R2*n
    if center:
        x0 -= 0.5*l*n
        x1 -= 0.5*l*n
        x2 -= 0.5*l*n
    l0 = Cylinder(scene, x0, x1, r)
    l1 = Cylinder(scene, x1 - 0.5*r*n, x1 - 0.5*r*n - R2*n + R1*t1, r)
    l2 = Cylinder(scene, x1 - 0.5*r*n, x1 - 0.5*r*n - R2*n - R1*t1, r)

    # Extract model
    return scene.to_model()

def UnitTetrahedron(color=(0.0, 1.0, 0.0, 0.5)):
    "Return model for unit tetrahedron."

    info("Plotting unit tetrahedron")

    # Create separate scene (since we will extract a model, not render)
    scene = soya.World()

    # Create vertices
    v0 = soya.Vertex(scene, 0.0, 0.0, 0.0, diffuse=color)
    v1 = soya.Vertex(scene, 1.0, 0.0, 0.0, diffuse=color)
    v2 = soya.Vertex(scene, 0.0, 1.0, 0.0, diffuse=color)
    v3 = soya.Vertex(scene, 0.0, 0.0, 1.0, diffuse=color)

    # Create edges
    e0 = Cylinder(scene, v0, v1, 0.005)
    e1 = Cylinder(scene, v0, v2, 0.005)
    e2 = Cylinder(scene, v0, v3, 0.005)
    e3 = Cylinder(scene, v1, v2, 0.005)
    e4 = Cylinder(scene, v1, v3, 0.005)
    e5 = Cylinder(scene, v2, v3, 0.005)

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

    info("Plotting unit triangle")

    # Create separate scene (since we will extract a model, not render)
    scene = soya.World()

    # Create vertice
    v0 = soya.Vertex(scene, 0.0, 0.0, 0.0, diffuse=color)
    v1 = soya.Vertex(scene, 1.0, 0.0, 0.0, diffuse=color)
    v2 = soya.Vertex(scene, 0.0, 1.0, 0.0, diffuse=color)

    # Create edges
    e0 = Cylinder(scene, v0, v1, 0.005)
    e1 = Cylinder(scene, v0, v2, 0.005)
    e2 = Cylinder(scene, v1, v2, 0.005)

    # Create face
    f = soya.Face(scene, (v0, v1, v2))

    # Make face double sided
    f.double_sided = 1

    # Extract model
    model = scene.to_model()

    return model

def PointEvaluation(x):
    "Return model for point evaluation at given point."

    info("Plotting dof: point evaluation at x = %s" % str(x))

    # Make sure point is 3D
    x = to3d(x)

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

def PointDerivative(x):
    "Return model for evaluation of derivatives at given point."

    info("Plotting dof: point derivative at x = %s" % str(x))

    # Make sure point is 3D
    x = to3d(x)

    # Create separate scene (since we will extract a model, not render)
    scene = soya.World()

    # Define material (color) for the sphere
    material = soya.Material()
    material.diffuse = (0.0, 0.0, 0.0, 0.2)

    # Create sphere
    sphere = Sphere(scene, material=material)

    # Scale and moveand move to coordinate
    sphere.scale(0.1, 0.1, 0.1)
    p = sphere.position()
    p.set_xyz(x[0], x[1], x[2])
    sphere.move(p)

    # Extract model
    model = scene.to_model()

    return model

def PointSecondDerivative(x):
    "Return model for evaluation of second derivatives at given point."

    info("Plotting dof: point derivative at x = %s" % str(x))

    # Make sure point is 3D
    x = to3d(x)

    # Create separate scene (since we will extract a model, not render)
    scene = soya.World()

    # Define material (color) for the sphere
    material = soya.Material()
    material.diffuse = (0.0, 0.0, 0.0, 0.05)

    # Create sphere
    sphere = Sphere(scene, material=material)

    # Scale and moveand move to coordinate
    sphere.scale(0.15, 0.15, 0.15)
    p = sphere.position()
    p.set_xyz(x[0], x[1], x[2])
    sphere.move(p)

    # Extract model
    model = scene.to_model()

    return model

def DirectionalEvaluation(x, n, flip=False, center=False):
    "Return model for directional evaluation at given point in given direction."

    info("Plotting dof: directional evaluation at x = %s in direction n = %s" % (str(x), str(n)))

    # Make sure points are 3D
    x = to3d(x)
    n = to3d(n)

    # Create separate scene (since we will extract a model, not render)
    scene = soya.World()

    # Normalize
    n = array(n)
    n = 0.75 * n / norm(n)

    # Flip normal if necessary
    if flip and not pointing_outwards(x, n):
        info("Flipping direction of arrow so it points outward.")
        n = -n

    # Create arrow
    arrow = Arrow(scene, x, n, center)

    # Extract model
    model = scene.to_model()

    return model

def DirectionalDerivative(x, n):
    "Return model for directional derivative at given point in given direction."

    info("Plotting dof: directional derivative at x = %s in direction n = %s" % (str(x), str(n)))

    # Make sure points are 3D
    x = to3d(x)
    n = to3d(n)

    # Create separate scene (since we will extract a model, not render)
    scene = soya.World()

    # Normalize
    n = array(n)
    n = 0.75 * n / norm(n)

    # Create line
    line = Cylinder(scene, x - 0.07*n, x + 0.07*n, 0.005)

    # Extract model
    model = scene.to_model()

    return model

def IntegralMoment(element, num_moments):
    "Return model for integral moment for given element."

    info("Plotting dof: integral moment")

    # Set position
    if element.cell().domain() == "triangle":
        a = 1.0 / (2 + sqrt(2)) # this was a fun exercise
        x = (a, a, 0.0)
    else:
        a = 1.0 / (3 + sqrt(3)) # so was this
        x = (a, a, a)

    # Fancy scaling of radius and color
    r = 1.0 / (num_moments + 5)
    if num_moments % 2 == 0:
        c = 1.0
    else:
        c = 0.0

    # Create separate scene (since we will extract a model, not render)
    scene = soya.World()

    # Define material (color) for the sphere
    material = soya.Material()
    material.diffuse = (c, c, c, 0.7)

    # Create sphere
    sphere = Sphere(scene, material=material)

    # Scale and moveand move to coordinate
    sphere.scale(r, r, r)
    p = sphere.position()
    p.set_xyz(x[0], x[1], x[2])
    sphere.move(p)

    # Extract model
    model = scene.to_model()

    return model

def create_cell_model(element):
    "Create Soya3D model for cell."

    # Get color
    family = element.family()
    if not family in element_colors:
        warning("Don't know a good color for elements of type '%s', using default color." % family)
        family = "Lagrange"
    color = element_colors[family]
    color = (color[0], color[1], color[2], 0.7)

    # Create model based on domain type
    domain = element.cell().domain()
    if domain == "triangle":
        return UnitTriangle(color), False
    elif domain == "tetrahedron":
        return UnitTetrahedron(color), True

    error("Unable to plot element, unhandled cell type: %s" % str(domain))

def create_dof_models(element):
    "Create Soya3D models for dofs."

    # Flags for whether to flip and center arrows
    directional = {"PointScaledNormalEval": (True,  False),
                   "PointEdgeTangent":      (False, True),
                   "PointFaceTangent":      (False, True)}

    # Elements not supported fully by FIAT
    unsupported = {"Argyris":            argyris_dofs,
                   "Arnold-Winther":     arnold_winther_dofs,
                   "Hermite":            hermite_dofs,
                   "Mardal-Tai-Winther": mardal_tai_winther_dofs,
                   "Morley":             morley_dofs}

    # Check if element is supported
    family = element.family()
    if not family in unsupported:
        # Create FIAT element and get dofs
        fiat_element = create_element(element)
        dofs = [(dof.get_type_tag(), dof.get_point_dict()) for dof in fiat_element.dual_basis()]

    else:

        # Bybass FIAT and set the dofs ourselves
        dofs = unsupported[family](element)

    # Iterate over dofs and add models
    models = []
    num_moments = 0
    for (dof_type, L) in dofs:

        # Check type of dof
        if dof_type == "PointEval":

            # Point evaluation, just get point
            points = L.keys()
            if not len(points) == 1:
                error("Strange dof, single point expected.")
            x = points[0]

            # Generate model
            models.append(PointEvaluation(x))

        elif dof_type == "PointNormalDeriv":

            # Evaluation of derivatives at point
            points = L.keys()
            if not len(points) == 1:
                error("Strange dof, single point expected.")
            x = points[0]
            n = [xx[0] for xx in L[x]]

            # Generate model
            models.append(DirectionalDerivative(x, n))

        elif dof_type == "PointDeriv":

            # Evaluation of derivatives at point
            points = L.keys()
            if not len(points) == 1:
                error("Strange dof, single point expected.")
            x = points[0]

            # Generate model
            models.append(PointDerivative(x))

        elif dof_type == "PointSecondDeriv":

            # Evaluation of derivatives at point
            points = L.keys()
            if not len(points) == 1:
                error("Strange dof, single point expected.")
            x = points[0]

            # Generate model
            models.append(PointSecondDerivative(x))

        elif dof_type in directional:

            # Normal evaluation, get point and normal
            points = L.keys()
            if not len(points) == 1:
                error("Strange dof, single point expected.")
            x = points[0]
            n = [xx[0] for xx in L[x]]

            # Generate model
            flip, center = directional[dof_type]
            models.append(DirectionalEvaluation(x, n, flip, center))

        elif dof_type in ("FrobeniusIntegralMoment", "IntegralMoment", "ComponentPointEval"):

            # Generate model
            models.append(IntegralMoment(element, num_moments))

            # Count the number of integral moments
            num_moments += 1

        else:
            error("Unable to plot dof, unhandled dof type: %s" % str(dof_type))

    return models, num_moments

def pointing_outwards(x, n):
    "Check if n is pointing inwards, used for flipping dofs."
    eps = 1e-10
    x = array(x) + 0.1*array(n)
    return x[0] < -eps or x[1] < -eps or x[2] < -eps or x[2] > 1.0 - x[0] - x[1] + eps

def to3d(x):
    "Make sure point is 3D."
    if len(x) == 2:
        x = (x[0], x[1], 0.0)
    return x

def arnold_winther_dofs(element):
    "Special fix for Arnold-Winther elements until Rob fixes in FIAT."

    if not element.cell().domain() == "triangle":
        error("Unable to plot element, only know how to plot Mardal-Tai-Winther on triangles.")

    return [("PointEval",        {(0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointEval",        {(0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointEval",        {(0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointEval",        {(1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointEval",        {(1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointEval",        {(1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointEval",        {(0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointEval",        {(0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointEval",        {(0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointScaledNormalEval", {(1.0/7, 0.0):     [  (0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(2.0/7, 0.0):     [  (0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(3.0/7, 0.0):     [  (0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(4.0/7, 0.0):     [  (0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(5.0/7, 0.0):     [  (0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(6.0/7, 0.0):     [  (0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(6.0/7, 1.0/7.0): [  (1.0, (0,)),  (1.0, (1,))]}),
            ("PointScaledNormalEval", {(5.0/7, 2.0/7.0): [  (1.0, (0,)),  (1.0, (1,))]}),
            ("PointScaledNormalEval", {(4.0/7, 3.0/7.0): [  (1.0, (0,)),  (1.0, (1,))]}),
            ("PointScaledNormalEval", {(3.0/7, 4.0/7.0): [  (1.0, (0,)),  (1.0, (1,))]}),
            ("PointScaledNormalEval", {(2.0/7, 5.0/7.0): [  (1.0, (0,)),  (1.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0,   1.0/7.0): [ (-1.0, (0,)),  (0.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0,   2.0/7.0): [ (-1.0, (0,)),  (0.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0,   3.0/7.0): [ (-1.0, (0,)),  (0.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0,   4.0/7.0): [ (-1.0, (0,)),  (0.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0,   5.0/7.0): [ (-1.0, (0,)),  (0.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0,   6.0/7.0): [ (-1.0, (0,)),  (0.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0,   6.0/7.0): [ (-1.0, (0,)),  (0.0, (1,))]}),
            ("IntegralMoment", None),
            ("IntegralMoment", None),
            ("IntegralMoment", None)]

def argyris_dofs(element):
    "Special fix for Hermite elements until Rob fixes in FIAT."

    if not element.degree() == 5:
        error("Unable to plot element, only know how to plot quintic Argyris elements.")

    if not element.cell().domain() == "triangle":
        error("Unable to plot element, only know how to plot Argyris on triangles.")

    return [("PointEval",        {(0.0, 0.0): [ (1.0, ()) ]}),
            ("PointEval",        {(1.0, 0.0): [ (1.0, ()) ]}),
            ("PointEval",        {(0.0, 1.0): [ (1.0, ()) ]}),
            ("PointDeriv",       {(0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof twice
            ("PointDeriv",       {(0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof twice
            ("PointDeriv",       {(1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof twice
            ("PointDeriv",       {(1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof twice
            ("PointDeriv",       {(0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof twice
            ("PointDeriv",       {(0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof twice
            ("PointSecondDeriv", {(0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointSecondDeriv", {(0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointSecondDeriv", {(0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointSecondDeriv", {(1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointSecondDeriv", {(1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointSecondDeriv", {(1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointSecondDeriv", {(0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointSecondDeriv", {(0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointSecondDeriv", {(0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof three times
            ("PointNormalDeriv", {(0.5, 0.0): [ (0.0, (0,)), (-1.0,  (1,))]}),
            ("PointNormalDeriv", {(0.5, 0.5): [ (1.0, (0,)), ( 1.0,  (1,))]}),
            ("PointNormalDeriv", {(0.0, 0.5): [(-1.0, (0,)), ( 0.0,  (1,))]})]

def hermite_dofs(element):
    "Special fix for Hermite elements until Rob fixes in FIAT."

    dofs_2d = [("PointEval",  {(0.0, 0.0): [ (1.0, ()) ]}),
               ("PointEval",  {(1.0, 0.0): [ (1.0, ()) ]}),
               ("PointEval",  {(0.0, 1.0): [ (1.0, ()) ]}),
               ("PointDeriv", {(0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof twice
               ("PointDeriv", {(0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof twice
               ("PointDeriv", {(1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof twice
               ("PointDeriv", {(1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof twice
               ("PointDeriv", {(0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof twice
               ("PointDeriv", {(0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof twice
               ("PointEval",  {(1.0/3, 1.0/3): [ (1.0, ()) ]})]

    dofs_3d = [("PointEval",  {(0.0, 0.0, 0.0): [ (1.0, ()) ]}),
               ("PointEval",  {(1.0, 0.0, 0.0): [ (1.0, ()) ]}),
               ("PointEval",  {(0.0, 1.0, 0.0): [ (1.0, ()) ]}),
               ("PointEval",  {(0.0, 0.0, 1.0): [ (1.0, ()) ]}),
               ("PointDeriv", {(0.0, 0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
               ("PointDeriv", {(0.0, 0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
               ("PointDeriv", {(0.0, 0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
               ("PointDeriv", {(1.0, 0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
               ("PointDeriv", {(1.0, 0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
               ("PointDeriv", {(1.0, 0.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
               ("PointDeriv", {(0.0, 1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
               ("PointDeriv", {(0.0, 1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
               ("PointDeriv", {(0.0, 1.0, 0.0): [ (1.0, ()) ]}), # hack, same dof three times
               ("PointDeriv", {(0.0, 0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof three times
               ("PointDeriv", {(0.0, 0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof three times
               ("PointDeriv", {(0.0, 0.0, 1.0): [ (1.0, ()) ]}), # hack, same dof three times
               ("PointEval",  {(1.0/3, 1.0/3, 1.0/3): [ (1.0, ()) ]}),
               ("PointEval",  {(0.0,   1.0/3, 1.0/3): [ (1.0, ()) ]}),
               ("PointEval",  {(1.0/3, 0.0,   1.0/3): [ (1.0, ()) ]}),
               ("PointEval",  {(1.0/3, 1.0/3, 0.0):   [ (1.0, ()) ]})]

    if element.cell().domain() == "triangle":
        return dofs_2d
    else:
        return dofs_3d

def mardal_tai_winther_dofs(element):
    "Special fix for Mardal-Tai-Winther elements until Rob fixes in FIAT."

    if not element.cell().domain() == "triangle":
        error("Unable to plot element, only know how to plot Mardal-Tai-Winther on triangles.")

    return [("PointScaledNormalEval", {(1.0/3, 0.0):     [  (0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(2.0/3, 0.0):     [  (0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(2.0/3, 1.0/3.0): [  (1.0, (0,)),  (1.0, (1,))]}),
            ("PointScaledNormalEval", {(1.0/3, 2.0/3.0): [  (1.0, (0,)),  (1.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0,   1.0/3.0): [ (-1.0, (0,)),  (0.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0,   2.0/3.0): [ (-1.0, (0,)),  (0.0, (1,))]}),
            ("PointEdgeTangent", {(0.5, 0.0): [ (-1.0, (0,)),  (0.0, (1,))]}),
            ("PointEdgeTangent", {(0.5, 0.5): [ (-1.0, (0,)),  (1.0, (1,))]}),
            ("PointEdgeTangent", {(0.0, 0.5): [  (0.0, (0,)), (-1.0, (1,))]})]

def morley_dofs(element):
    "Special fix for Morley elements until Rob fixes in FIAT."

    if not element.cell().domain() == "triangle":
        error("Unable to plot element, only know how to plot Morley on triangles.")

    return [("PointEval",        {(0.0, 0.0): [ (1.0, ()) ]}),
            ("PointEval",        {(1.0, 0.0): [ (1.0, ()) ]}),
            ("PointEval",        {(0.0, 1.0): [ (1.0, ()) ]}),
            ("PointNormalDeriv", {(0.5, 0.0): [ (0.0, (0,)), (-1.0,  (1,))]}),
            ("PointNormalDeriv", {(0.5, 0.5): [ (1.0, (0,)), ( 1.0,  (1,))]}),
            ("PointNormalDeriv", {(0.0, 0.5): [(-1.0, (0,)), ( 0.0,  (1,))]})]
