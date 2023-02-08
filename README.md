# timex_hackathon
Toy model for timex hackathon.

Contact: Erik Schnaubelt, erik.schnaubelt@cern.ch

## Geometry/Mesh

Main file is `pancake_builder.py`. It uses `inputs/thermal_coil_TSA.yaml`.
Calling `pancake_builder` creates the geometry and mesh of the model in the
folder `thermal_coil_TSA`.

## FE Solution

Main file is `thermal.pro`. Can be run via CLI or GUI.

## Requirements: Gmsh/Getdp

In order to run the FE simulation, executables of both Gmsh and GetDP are required.
For Gmsh, please download it from [the Gmsh website](https://gmsh.info/).
For GetDP, a specific compiled version of GetDP is required to have access to
all material functions. Please download it from the
[the CERN Gitlab](https://gitlab.cern.ch/steam/cerngetdp/-/releases/cerngetdp_v1.05)
under `Assets/Other`.
Please note that at minimum CERNGetDP 1.05 is required.
