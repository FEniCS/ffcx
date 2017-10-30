# -*- coding: utf-8 -*-
"This module provides functionality for plotting finite elements."

# Copyright (C) 2010 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2010-12-07
# Last changed: 2010-12-15

__all__ = ["plot"]

from numpy import dot, cross, array, sin, cos, pi, sqrt
from numpy.linalg import norm
import sys

from ffc.fiatinterface import create_element
from ffc.log import warning, error, info

# Import Soya3D
try:
    import soya
    from soya.sphere import Sphere
    from soya.label3d import Label3D
    from soya.sdlconst import QUIT
    _soya_imported = True
except ImportError:
    _soya_imported = False

# Colors for elements
element_colors = {"Argyris": (0.45, 0.70, 0.80),
                  "Arnold-Winther": (0.00, 0.00, 1.00),
                  "Brezzi-Douglas-Marini": (1.00, 1.00, 0.00),
                  "Crouzeix-Raviart": (1.00, 0.25, 0.25),
                  "Discontinuous Lagrange": (0.00, 0.25, 0.00),
                  "Discontinuous Raviart-Thomas": (0.90, 0.90, 0.30),
                  "Hermite": (0.50, 1.00, 0.50),
                  "Lagrange": (0.00, 1.00, 0.00),
                  "Mardal-Tai-Winther": (1.00, 0.10, 0.90),
                  "Morley": (0.40, 0.40, 0.40),
                  "Nedelec 1st kind H(curl)": (0.90, 0.30, 0.00),
                  "Nedelec 2nd kind H(curl)": (0.70, 0.20, 0.00),
                  "Raviart-Thomas": (0.90, 0.60, 0.00)}


def plot(element, rotate=True):
    "Plot finite element."

    # Check if Soya3D has been imported
    if not _soya_imported:
        warning("Unable to plot element, Soya3D not available (install package python-soya).")
        return

    # Special case: plot dof notation
    if element == "notation":

        # Create model for notation
        notation = create_notation_models()

        # Render plot window
        render(notation, "Notation", 0, True, rotate)

    else:

        # Create cell model
        cell, is3d = create_cell_model(element)

        cellname = element.cell().cellname()  # Assuming single cell

        # Create dof models
        dofs, num_moments = create_dof_models(element)

        # Create title
        if element.degree() is not None:
            title = "%s of degree %d on a %s" % (element.family(), element.degree(), cellname)
        else:
            title = "%s on a %s" % (element.family(), cellname)

        # Render plot window
        render([cell] + dofs, title, num_moments, is3d, rotate)


def render(models, title, num_moments, is3d, rotate):
    "Render given list of models."

    # Note that we view from the positive z-axis, and not from the
    # negative y-axis. This should make no difference since the
    # element dofs are symmetric anyway and it plays better with
    # the default camera settings in Soya.

    # Initialize Soya
    soya.init(title)

    # Create scene
    scene = soya.World()
    scene.atmosphere = soya.Atmosphere()
    if title == "Notation":
        scene.atmosphere.bg_color = (0.0, 1.0, 0.0, 1.0)
    else:
        scene.atmosphere.bg_color = (1.0, 1.0, 1.0, 1.0)

    # Not used, need to manually handle rotation
    # label = Label3D(scene, text=str(num_moments), size=0.005)
    # label.set_xyz(1.0, 1.0, 1.0)
    # label.set_color((0.0, 0.0, 0.0, 1.0))

    # Define rotation
    if is3d:
        class RotatingBody(soya.Body):

            def advance_time(self, proportion):
                self.rotate_y(2.0 * proportion)
    else:
        class RotatingBody(soya.Body):

            def advance_time(self, proportion):
                self.rotate_z(2.0 * proportion)

    # Select type of display, rotating or not
    if rotate:
        Body = RotatingBody
    else:
        Body = soya.Body

    # Add all models
    for model in models:
        body = Body(scene, model)

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
    p = camera.position()
    if is3d:
        if rotate:
            camera.set_xyz(-20, 10, 50.0)
            camera.fov = 2.1
            p.set_xyz(0.0, 0.4, 0.0)
        else:
            camera.set_xyz(-20, 10, 50.0)
            camera.fov = 1.6
            p.set_xyz(0.3, 0.42, 0.5)
    else:
        if rotate:
            camera.set_xyz(0, 10, 50.0)
            camera.fov = 2.6
            p.set_xyz(0.0, 0.0, 0.0)
        else:
            camera.set_xyz(0, 10, 50.0)
            camera.fov = 1.7
            p.set_xyz(0.5, 0.4, 0.0)
    camera.look_at(p)
    soya.set_root_widget(camera)

    # Handle exit
    class Idler(soya.Idler):

        def end_round(self):
            for event in self.events:
                if event[0] == QUIT:
                    print("Closing plot, bye bye")
                    sys.exit(0)

    # Main loop
    idler = Idler(scene)
    idler.idle()


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
    dtheta = 2.0 * pi / float(num_steps)
    for i in range(num_steps):

        # Compute coordinates for square
        dx0 = cos(i * dtheta) * t0 + sin(i * dtheta) * t1
        dx1 = cos((i + 1) * dtheta) * t0 + sin((i + 1) * dtheta) * t1
        x0 = p0 + r * dx0
        x1 = p0 + r * dx1
        x2 = p1 + r * dx0
        x3 = p1 + r * dx1

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


def Cone(scene, p0, p1, r, color=(0.0, 0.0, 0.0, 1.0)):
    "Return model for cone from p0 to p1 with radius r."

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
    dtheta = 2.0 * pi / float(num_steps)
    v2 = soya.Vertex(scene, p1[0], p1[1], p1[2], diffuse=color)
    for i in range(num_steps):

        # Compute coordinates for bottom of face
        dx0 = cos(i * dtheta) * t0 + sin(i * dtheta) * t1
        dx1 = cos((i + 1) * dtheta) * t0 + sin((i + 1) * dtheta) * t1
        x0 = p0 + r * dx0
        x1 = p0 + r * dx1

        # Create face
        v0 = soya.Vertex(scene, x0[0], x0[1], x0[2], diffuse=color)
        v1 = soya.Vertex(scene, x1[0], x1[1], x1[2], diffuse=color)
        f = soya.Face(scene, (v0, v1, v2))
        f.double_sided = 1

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
    L = 0.3
    l = 0.35 * L  # noqa: E741
    r = 0.04 * L
    R = 0.125 * L

    # Center arrow
    if center:
        print("Centering!")
        x -= 0.5 * (L + l) * n

    # Create cylinder and cone
    cylinder = Cylinder(scene, x, x + L * n, r)
    cone = Cone(scene, x + L * n, x + (L + l) * n, R)

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
    e0 = Cylinder(scene, v0, v1, 0.007)
    e1 = Cylinder(scene, v0, v2, 0.007)
    e2 = Cylinder(scene, v0, v3, 0.007)
    e3 = Cylinder(scene, v1, v2, 0.007)
    e4 = Cylinder(scene, v1, v3, 0.007)
    e5 = Cylinder(scene, v2, v3, 0.007)

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
    e0 = Cylinder(scene, v0, v1, 0.007)
    e1 = Cylinder(scene, v0, v2, 0.007)
    e2 = Cylinder(scene, v1, v2, 0.007)

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
    line = Cylinder(scene, x - 0.07 * n, x + 0.07 * n, 0.005)

    # Extract model
    model = scene.to_model()

    return model


def IntegralMoment(cellname, num_moments, x=None):
    "Return model for integral moment for given element."

    info("Plotting dof: integral moment")

    # Set position
    if x is None and cellname == "triangle":
        a = 1.0 / (2 + sqrt(2))  # this was a fun exercise
        x = (a, a, 0.0)
    elif x is None:
        a = 1.0 / (3 + sqrt(3))  # so was this
        x = (a, a, a)

    # Make sure point is 3D
    x = to3d(x)

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
    if family not in element_colors:
        warning("Don't know a good color for elements of type '%s', using default color." % family)
        family = "Lagrange"
    color = element_colors[family]
    color = (color[0], color[1], color[2], 0.7)

    # Create model based on cell type
    cellname = element.cell().cellname()
    if cellname == "triangle":
        return UnitTriangle(color), False
    elif cellname == "tetrahedron":
        return UnitTetrahedron(color), True

    error("Unable to plot element, unhandled cell type: %s" % str(cellname))


def create_dof_models(element):
    "Create Soya3D models for dofs."

    # Flags for whether to flip and center arrows
    directional = {"PointScaledNormalEval": (True, False),
                   "PointEdgeTangent": (False, True),
                   "PointFaceTangent": (False, True)}

    # Elements not supported fully by FIAT
    unsupported = {"Argyris": argyris_dofs,
                   "Arnold-Winther": arnold_winther_dofs,
                   "Hermite": hermite_dofs,
                   "Mardal-Tai-Winther": mardal_tai_winther_dofs,
                   "Morley": morley_dofs}

    # Check if element is supported
    family = element.family()
    if family not in unsupported:
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
            points = list(L.keys())
            if not len(points) == 1:
                error("Strange dof, single point expected.")
            x = points[0]

            # Generate model
            models.append(PointEvaluation(x))

        elif dof_type == "PointDeriv":

            # Evaluation of derivatives at point
            points = list(L.keys())
            if not len(points) == 1:
                error("Strange dof, single point expected.")
            x = points[0]

            # Generate model
            models.append(PointDerivative(x))

        elif dof_type == "PointSecondDeriv":

            # Evaluation of derivatives at point
            points = list(L.keys())
            if not len(points) == 1:
                error("Strange dof, single point expected.")
            x = points[0]

            # Generate model
            models.append(PointSecondDerivative(x))

        elif dof_type in directional:

            # Normal evaluation, get point and normal
            points = list(L.keys())
            if not len(points) == 1:
                error("Strange dof, single point expected.")
            x = points[0]
            n = [xx[0] for xx in L[x]]

            # Generate model
            flip, center = directional[dof_type]
            models.append(DirectionalEvaluation(x, n, flip, center))

        elif dof_type == "PointNormalDeriv":

            # Evaluation of derivatives at point
            points = list(L.keys())
            if not len(points) == 1:
                error("Strange dof, single point expected.")
            x = points[0]
            n = [xx[0] for xx in L[x]]

            # Generate model
            models.append(DirectionalDerivative(x, n))

        elif dof_type in ("FrobeniusIntegralMoment", "IntegralMoment", "ComponentPointEval"):

            # Generate model
            models.append(IntegralMoment(element.cell().cellname(), num_moments))

            # Count the number of integral moments
            num_moments += 1

        else:
            error("Unable to plot dof, unhandled dof type: %s" % str(dof_type))

    return models, num_moments


def create_notation_models():
    "Create Soya 3D models for notation."

    models = []

    y = 1.3
    dy = -0.325

    # Create model for evaluation
    models.append(PointEvaluation([0, y]))
    y += dy

    # Create model for derivative evaluation
    models.append(PointDerivative([0, y]))
    models.append(PointDerivative([0, y]))
    models.append(PointDerivative([0, y]))
    y += dy

    # Create model for second derivative evaluation
    models.append(PointSecondDerivative([0, y]))
    models.append(PointSecondDerivative([0, y]))
    models.append(PointSecondDerivative([0, y]))
    y += dy

    # Create model for directional evaluation
    models.append(DirectionalEvaluation([0, y], [1, 1], False, True))
    y += dy

    # Create model for directional evaluation
    models.append(DirectionalDerivative([0, y], [1, 1]))
    y += dy

    # Create model for integral moments
    models.append(IntegralMoment("tetrahedron", 0, [0, y]))
    models.append(IntegralMoment("tetrahedron", 1, [0, y]))
    models.append(IntegralMoment("tetrahedron", 2, [0, y]))

    return models


def pointing_outwards(x, n):
    "Check if n is pointing inwards, used for flipping dofs."
    eps = 1e-10
    x = array(x) + 0.1 * array(n)
    return x[0] < -eps or x[1] < -eps or x[2] < -eps or x[2] > 1.0 - x[0] - x[1] + eps


def to3d(x):
    "Make sure point is 3D."
    if len(x) == 2:
        x = (x[0], x[1], 0.0)
    return x


def arnold_winther_dofs(element):
    "Special fix for Arnold-Winther elements until Rob fixes in FIAT."

    if not element.cell().cellname() == "triangle":
        error("Unable to plot element, only know how to plot Mardal-Tai-Winther on triangles.")

    return [("PointEval", {(0.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointEval", {(0.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointEval", {(0.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointEval", {(1.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointEval", {(1.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointEval", {(1.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointEval", {(0.0, 1.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointEval", {(0.0, 1.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointEval", {(0.0, 1.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointScaledNormalEval", {(1.0 / 5, 0.0): [(0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(2.0 / 5, 0.0): [(0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(3.0 / 5, 0.0): [(0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(4.0 / 5, 0.0): [(0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(4.0 / 5, 1.0 / 5.0): [(1.0, (0,)), (1.0, (1,))]}),
            ("PointScaledNormalEval", {(3.0 / 5, 2.0 / 5.0): [(1.0, (0,)), (1.0, (1,))]}),
            ("PointScaledNormalEval", {(2.0 / 5, 3.0 / 5.0): [(1.0, (0,)), (1.0, (1,))]}),
            ("PointScaledNormalEval", {(1.0 / 5, 4.0 / 5.0): [(1.0, (0,)), (1.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0, 1.0 / 5.0): [(-1.0, (0,)), (0.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0, 2.0 / 5.0): [(-1.0, (0,)), (0.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0, 3.0 / 5.0): [(-1.0, (0,)), (0.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0, 4.0 / 5.0): [(-1.0, (0,)), (0.0, (1,))]}),
            ("IntegralMoment", None),
            ("IntegralMoment", None),
            ("IntegralMoment", None)]


def argyris_dofs(element):
    "Special fix for Hermite elements until Rob fixes in FIAT."

    if not element.degree() == 5:
        error("Unable to plot element, only know how to plot quintic Argyris elements.")

    if not element.cell().cellname() == "triangle":
        error("Unable to plot element, only know how to plot Argyris on triangles.")

    return [("PointEval", {(0.0, 0.0): [(1.0, ())]}),
            ("PointEval", {(1.0, 0.0): [(1.0, ())]}),
            ("PointEval", {(0.0, 1.0): [(1.0, ())]}),
            ("PointDeriv", {(0.0, 0.0): [(1.0, ())]}),  # hack, same dof twice
            ("PointDeriv", {(0.0, 0.0): [(1.0, ())]}),  # hack, same dof twice
            ("PointDeriv", {(1.0, 0.0): [(1.0, ())]}),  # hack, same dof twice
            ("PointDeriv", {(1.0, 0.0): [(1.0, ())]}),  # hack, same dof twice
            ("PointDeriv", {(0.0, 1.0): [(1.0, ())]}),  # hack, same dof twice
            ("PointDeriv", {(0.0, 1.0): [(1.0, ())]}),  # hack, same dof twice
            ("PointSecondDeriv", {(0.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointSecondDeriv", {(0.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointSecondDeriv", {(0.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointSecondDeriv", {(1.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointSecondDeriv", {(1.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointSecondDeriv", {(1.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointSecondDeriv", {(0.0, 1.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointSecondDeriv", {(0.0, 1.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointSecondDeriv", {(0.0, 1.0): [(1.0, ())]}),  # hack, same dof three times
            ("PointNormalDeriv", {(0.5, 0.0): [(0.0, (0,)), (-1.0, (1,))]}),
            ("PointNormalDeriv", {(0.5, 0.5): [(1.0, (0,)), (1.0, (1,))]}),
            ("PointNormalDeriv", {(0.0, 0.5): [(-1.0, (0,)), (0.0, (1,))]})]


def hermite_dofs(element):
    "Special fix for Hermite elements until Rob fixes in FIAT."

    dofs_2d = [("PointEval", {(0.0, 0.0): [(1.0, ())]}),
               ("PointEval", {(1.0, 0.0): [(1.0, ())]}),
               ("PointEval", {(0.0, 1.0): [(1.0, ())]}),
               ("PointDeriv", {(0.0, 0.0): [(1.0, ())]}),  # hack, same dof twice
               ("PointDeriv", {(0.0, 0.0): [(1.0, ())]}),  # hack, same dof twice
               ("PointDeriv", {(1.0, 0.0): [(1.0, ())]}),  # hack, same dof twice
               ("PointDeriv", {(1.0, 0.0): [(1.0, ())]}),  # hack, same dof twice
               ("PointDeriv", {(0.0, 1.0): [(1.0, ())]}),  # hack, same dof twice
               ("PointDeriv", {(0.0, 1.0): [(1.0, ())]}),  # hack, same dof twice
               ("PointEval", {(1.0 / 3, 1.0 / 3): [(1.0, ())]})]

    dofs_3d = [("PointEval", {(0.0, 0.0, 0.0): [(1.0, ())]}),
               ("PointEval", {(1.0, 0.0, 0.0): [(1.0, ())]}),
               ("PointEval", {(0.0, 1.0, 0.0): [(1.0, ())]}),
               ("PointEval", {(0.0, 0.0, 1.0): [(1.0, ())]}),
               ("PointDeriv", {(0.0, 0.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
               ("PointDeriv", {(0.0, 0.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
               ("PointDeriv", {(0.0, 0.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
               ("PointDeriv", {(1.0, 0.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
               ("PointDeriv", {(1.0, 0.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
               ("PointDeriv", {(1.0, 0.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
               ("PointDeriv", {(0.0, 1.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
               ("PointDeriv", {(0.0, 1.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
               ("PointDeriv", {(0.0, 1.0, 0.0): [(1.0, ())]}),  # hack, same dof three times
               ("PointDeriv", {(0.0, 0.0, 1.0): [(1.0, ())]}),  # hack, same dof three times
               ("PointDeriv", {(0.0, 0.0, 1.0): [(1.0, ())]}),  # hack, same dof three times
               ("PointDeriv", {(0.0, 0.0, 1.0): [(1.0, ())]}),  # hack, same dof three times
               ("PointEval", {(1.0 / 3, 1.0 / 3, 1.0 / 3): [(1.0, ())]}),
               ("PointEval", {(0.0, 1.0 / 3, 1.0 / 3): [(1.0, ())]}),
               ("PointEval", {(1.0 / 3, 0.0, 1.0 / 3): [(1.0, ())]}),
               ("PointEval", {(1.0 / 3, 1.0 / 3, 0.0): [(1.0, ())]})]

    if element.cell().cellname() == "triangle":
        return dofs_2d
    else:
        return dofs_3d


def mardal_tai_winther_dofs(element):
    "Special fix for Mardal-Tai-Winther elements until Rob fixes in FIAT."

    if not element.cell().cellname() == "triangle":
        error("Unable to plot element, only know how to plot Mardal-Tai-Winther on triangles.")

    return [("PointScaledNormalEval", {(1.0 / 3, 0.0): [(0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(2.0 / 3, 0.0): [(0.0, (0,)), (-1.0, (1,))]}),
            ("PointScaledNormalEval", {(2.0 / 3, 1.0 / 3.0): [(1.0, (0,)), (1.0, (1,))]}),
            ("PointScaledNormalEval", {(1.0 / 3, 2.0 / 3.0): [(1.0, (0,)), (1.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0, 1.0 / 3.0): [(-1.0, (0,)), (0.0, (1,))]}),
            ("PointScaledNormalEval", {(0.0, 2.0 / 3.0): [(-1.0, (0,)), (0.0, (1,))]}),
            ("PointEdgeTangent", {(0.5, 0.0): [(-1.0, (0,)), (0.0, (1,))]}),
            ("PointEdgeTangent", {(0.5, 0.5): [(-1.0, (0,)), (1.0, (1,))]}),
            ("PointEdgeTangent", {(0.0, 0.5): [(0.0, (0,)), (-1.0, (1,))]})]


def morley_dofs(element):
    "Special fix for Morley elements until Rob fixes in FIAT."

    if not element.cell().cellname() == "triangle":
        error("Unable to plot element, only know how to plot Morley on triangles.")

    return [("PointEval", {(0.0, 0.0): [(1.0, ())]}),
            ("PointEval", {(1.0, 0.0): [(1.0, ())]}),
            ("PointEval", {(0.0, 1.0): [(1.0, ())]}),
            ("PointNormalDeriv", {(0.5, 0.0): [(0.0, (0,)), (-1.0, (1,))]}),
            ("PointNormalDeriv", {(0.5, 0.5): [(1.0, (0,)), (1.0, (1,))]}),
            ("PointNormalDeriv", {(0.0, 0.5): [(-1.0, (0,)), (0.0, (1,))]})]
