# NASA's Plot3D to Gmsh's MSH converter

Small Python script to convert from Plot3D into Gmsh mesh format. In general
you'll need mesh in Plot3D format (script was tested with files from
<http://turbmodels.larc.nasa.gov>) and Neutral Map File for boundary
description (also available from NASA for 3D versions of the meshes).

## Usage

```sh
$ p3d2gmsh.py [-o OUT_FILE] [-m MAP_FILE] P3D_FILE

    P3D_FILE: name of file with the mesh
    MAP_FILE: Neutral Map File (if omitted script will search for a file with
              nmf extension and the name of mesh file)
    OUT_FILE: name of output file (if omitted script will save a file with msh
              extension and the name of the mesh file)
```
