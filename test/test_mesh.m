% TEST_MESH: Test a single mesh and it's property

clearvars, close all
clc

% Add auxiliary function and meshes
addpath('../mesh/');
addpath('../pre-processing/');

% Read mesh
tic
mesh = read_gmsh('step_3d.msh');
fprintf('Time to read the mesh = %f\n', toc)

% Plot mesh
tic
plot_mesh(mesh);
fprintf('Time to read the mesh = %f\n', toc)