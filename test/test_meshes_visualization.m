% TEST_MESH_VISUALIZATION: Compare numerous meshes

clearvars, close all
clc

addpath('../mesh/');
addpath('../pre-processing/');

mesh = read_gmsh('circle_2d.msh'); plot_mesh(mesh);
mesh = read_gmsh('cylinder_2d.msh'); plot_mesh(mesh);
mesh = read_gmsh('cylinder_3d.msh'); plot_mesh(mesh);
mesh = read_gmsh('ell_2d.msh'); plot_mesh(mesh);
mesh = read_gmsh('rectangle.msh'); plot_mesh(mesh);
mesh = read_gmsh('step_2d.msh'); plot_mesh(mesh);
mesh = read_gmsh('step_3d.msh'); plot_mesh(mesh);
mesh = read_gmsh('t1.msh'); plot_mesh(mesh);
mesh = read_gmsh('t13_data.msh'); plot_mesh(mesh);
mesh = read_gmsh('tiny.msh'); plot_mesh(mesh);
mesh = read_gmsh('view5.msh'); plot_mesh(mesh);