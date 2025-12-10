# wing_geometry_python
Creation of the geometry of one wing starting from its NACA designation -> Naca 4 letters
This repository uses python and a gmsh environment
# Parameters to modify in independent wing : 
- WING_NACA -> NACA 4 letters
- WING_ROOT_CHORD -> root chord
- WING_TIP_CHORD -> tip chord
- WING_SEMI_SPAN -> semi span
- WING_NSPAN -> number of spanwise divisions for the wing
- WING_CHORD_DIV-> number of points along the chord
- WING_POS_X -> x position of the wing (when attach to the fuselage)
# Parameters to modify in wing+nacelle:
- WING_NACA -> NACA 4 letters
- WING_ROOT_CHORD -> root chord
- WING_TIP_CHORD -> tip chord
- WING_SEMI_SPAN -> semi span
- WING_NSPAN -> number of spanwise divisions for the wing
- WING_CHORD_DIV-> number of points along the chord
- WING_POS_X -> x position of the wing (when attach to the fuselage)
- PYLON_LENGTH -> lengt of the pylon (vertical in y direction)
- PYLON_RADIUS -> radius of the pylon
- NACELLE_RADIUS -> radius of the nacelle
- NACELLE_LENGTH -> length of the nacelle
- mid_span_z = semi_span * 0.5  -> modify the position of the nacelle along the wing (z directio)

# Python package to install:
- pip install gmsh
- pip install numpy




# Author 
Mihai Serbu
