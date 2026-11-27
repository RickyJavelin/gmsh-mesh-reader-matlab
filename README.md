# gmsh-mesh-reader-matlab

A lightweight MATLAB library to read, visualize, and analyze 2D and 3D mesh files (`.msh`) generated with [Gmsh](https://gmsh.info/).

## Features

* **Mesh Parsing:** Read node coordinates, element connectivity, and physical tags from `.msh` files (`read_gmsh.m`).
* **Visualization:** Plot 2D and 3D meshes (`plot_mesh.m`).
* **Geometric & Quality Analysis:** Compute element areas (`compute_area.m`) and evaluate mesh quality metrics (`mesh_quality.m`).
* **Sample Meshes Included:** Ready-to-use 2D and 3D sample `.msh` geometries for quick testing.

## Repository Structure

```text
├── mesh/                   # Sample 2D and 3D Gmsh files (.msh)
│   ├── circle_2d.msh
│   ├── cylinder_2d.msh
│   ├── cylinder_3d.msh
│   ├── ell_2d.msh
│   ├── rectangle.msh
│   ├── step_2d.msh
│   ├── step_3d.msh
│   └── ...
├── pre-processing/         # Core MATLAB functions
│   ├── compute_area.m      # Computes element areas
│   ├── mesh_quality.m      # Evaluates mesh quality metrics
│   ├── plot_mesh.m         # Visualizes 2D/3D meshes
│   └── read_gmsh.m         # Parses Gmsh .msh files
└── test/                   # Example scripts and test cases
    ├── test_mesh.m
    └── test_meshes_visualization.m
