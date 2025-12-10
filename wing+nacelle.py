#!/usr/bin/env python3
# Copyright (c) 2025 [Mihai]
#
# This code is licensed under the MIT License.
# See the LICENSE file in the project root for full license terms.
import gmsh
import numpy as np
from math import pi, sqrt, cos, sin, atan

# ------------------- NACA 4-digit airfoil -------------------

def parse_naca4(code: str):
    if len(code) != 4 or not code.isdigit():
        raise ValueError("Codice NACA non valido")
    m = int(code[0]) / 100.0
    p = int(code[1]) / 10.0
    t = int(code[2:]) / 100.0
    return m, p, t

def naca4(m, p, t, c=1.0, n_points=80):
    beta = np.linspace(0, pi, n_points)
    x = 0.5 * (1 - np.cos(beta)) * c
    yt = 5 * t * c * (
        0.2969 * np.sqrt(x / c)
        - 0.1260 * (x / c)
        - 0.3516 * (x / c) ** 2
        + 0.2843 * (x / c) ** 3
        - 0.1015 * (x / c) ** 4
    )
    yc = np.zeros_like(x)
    dyc = np.zeros_like(x)
    if p > 0:
        for i, xi in enumerate(x):
            xi_c = xi / c
            if xi_c < p:
                yc[i] = m / p**2 * (2 * p * xi_c - xi_c**2) * c
                dyc[i] = 2 * m / p**2 * (p - xi_c)
            else:
                yc[i] = m / (1 - p)**2 * ((1 - 2*p) + 2*p*xi_c - xi_c**2) * c
                dyc[i] = 2 * m / (1 - p)**2 * (p - xi_c)
    theta = np.arctan(dyc)
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)
    x_all = np.concatenate((xl[::-1], xu[1:]))
    y_all = np.concatenate((yl[::-1], yu[1:]))
    return list(zip(x_all, y_all))

# ------------------- Wing geometry -------------------

def build_structured_wing_plane_entities(params):
    naca_code = params.get("WING_NACA", "2412")
    root_chord = params.get("WING_ROOT_CHORD", 2.0)
    tip_chord = params.get("WING_TIP_CHORD", 0.8)
    semi_span = params.get("WING_SEMI_SPAN", 8.0)
    n_span_pts = params.get("WING_NSPAN", 12)
    n_chord_pts = params.get("WING_CHORD_DIV", 50)
    wing_pos_x = params.get("WING_POS_X", 0.0)

    m, p, t = parse_naca4(naca_code)
    span_positions = np.linspace(0.0, semi_span, n_span_pts)

    wires = []
    all_profiles_pts = []

    for z in span_positions:
        rel = z / semi_span if semi_span != 0 else 0.0
        chord = root_chord + (tip_chord - root_chord) * rel
        profile = naca4(m, p, t, c=chord, n_points=max(20, n_chord_pts))
        pts = []
        for x, y in profile:
            px = wing_pos_x + x
            py = y
            pz = z
            pts.append(gmsh.model.occ.addPoint(px, py, pz))
        pts.append(pts[0])  # close loop
        all_profiles_pts.append(pts)
        curve = gmsh.model.occ.addBSpline(pts)
        wire = gmsh.model.occ.addCurveLoop([curve])
        wires.append(wire)

    # surfaces between profiles
    surfaces = []
    for i in range(len(wires) - 1):
        s = gmsh.model.occ.addThruSections([wires[i], wires[i+1]], makeSolid=False, makeRuled=True)
        for dim, tag in s:
            if dim == 2:
                surfaces.append(tag)

    gmsh.model.occ.synchronize()

    # spanwise splines for mesh guidance
    n_nodes_per_profile = len(all_profiles_pts[0])
    spanwise_splines = []
    for i in range(n_nodes_per_profile):
        node_column = [profile[i] for profile in all_profiles_pts]
        spline = gmsh.model.occ.addBSpline(node_column)
        spanwise_splines.append(spline)

    gmsh.model.occ.synchronize()

    wing_entities = surfaces + spanwise_splines
    return wing_entities

# ------------------- Nacelle helper -------------------

def add_circle_ring(x_center, radius, n_points):
    pts = []
    for i in range(n_points):
        angle = 2*pi*i/n_points
        px = x_center
        py = radius*cos(angle)
        pz = radius*sin(angle)
        pts.append(gmsh.model.occ.addPoint(px, py, pz))
    return pts

# ------------------- Nacelles -------------------

def build_single_nacelle(params):
    surf_nacelle = []
    vol_pylon = []

    # Determine mid-span position
    semi_span = params.get("WING_SEMI_SPAN", 8.0)
    mid_span_z = semi_span * 0.5   #change here to move nacelle and pylon along the wing 
    root_chord = params.get("WING_ROOT_CHORD", 2.0)
    tip_chord = params.get("WING_TIP_CHORD", 0.8)

    # linear chord interpolation
    rel = mid_span_z / semi_span if semi_span != 0 else 0.0
    chord_mid = root_chord + (tip_chord - root_chord) * rel
    wing_pos_x = params.get("WING_POS_X", 0.0)
    nacelle_center_x = wing_pos_x + chord_mid * 0.5  # center of chord

    # create cylindrical nacelle centered at chord mid-point
    rings = [
        add_circle_ring(nacelle_center_x - 0.5*params["NACELLE_LENGTH"], params["NACELLE_RADIUS"], 12),
        add_circle_ring(nacelle_center_x + 0.5*params["NACELLE_LENGTH"], params["NACELLE_RADIUS"], 12)
    ]

    for i in range(12):
        a, b, c, d = rings[0][i], rings[0][(i + 1) % 12], rings[1][(i + 1) % 12], rings[1][i]
        l1 = gmsh.model.occ.addLine(a, b)
        l2 = gmsh.model.occ.addLine(b, c)
        l3 = gmsh.model.occ.addLine(c, d)
        l4 = gmsh.model.occ.addLine(d, a)
        loop = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
        s = gmsh.model.occ.addPlaneSurface([loop])
        surf_nacelle.append(s)

    gmsh.model.occ.synchronize()

    # vertical placement below wing
    vertical_offset = -params["PYLON_LENGTH"] - params["NACELLE_RADIUS"]
    gmsh.model.occ.translate([(2, s) for s in surf_nacelle], 0, vertical_offset, mid_span_z+params["NACELLE_RADIUS"])

    # create pylon connecting nacelle to wing plane
    top = (nacelle_center_x, -0.04, mid_span_z + params["NACELLE_RADIUS"]) # modify the y value here to ensure the top touches the wing
    attach = (nacelle_center_x, vertical_offset + params["NACELLE_RADIUS"], mid_span_z)  # wing plane
    vec = (0, attach[1]-top[1], 0)
    pylon = gmsh.model.occ.addCylinder(*top, *vec, params["PYLON_RADIUS"])
    vol_pylon.append(pylon)

    gmsh.model.occ.synchronize()
    return surf_nacelle + vol_pylon





# ------------------- Main -------------------

def main():
    params = {
        "WING_NACA": "2412",
        "WING_ROOT_CHORD": 2.0,
        "WING_TIP_CHORD": 0.8,
        "WING_SEMI_SPAN": 8.0,
        "WING_NSPAN": 21,
        "WING_CHORD_DIV": 41,
        "WING_POS_X": 0.0,
        "NACELLE_POS_X": 0.5,
        "NACELLE_RADIUS": 0.3,
        "NACELLE_LENGTH": 1.2,
        "PYLON_RADIUS": 0.05,
        "PYLON_LENGTH": 0.3
        
    }

    gmsh.initialize()
    gmsh.model.add("aircraft_geometry")

    wing_entities = build_structured_wing_plane_entities(params)
    nacelle_entities = build_single_nacelle(params)

    # Here you could generate meshes if desired
    # gmsh.model.mesh.generate(2)

    gmsh.fltk.run()
    gmsh.finalize()

if __name__ == "__main__":
    main()

