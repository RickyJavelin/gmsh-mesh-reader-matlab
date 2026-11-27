function plot_mesh(mesh)
% PLOT_MESH Automatically visualizes 1D, 2D, or 3D unstructured mesh data structures.
%
% Syntax:
%   plot_mesh(mesh)
%
% Description:
%   Renders computational meshes of arbitrary spatial dimensions (1D, 2D, or 3D) 
%   by parsing nodal coordinates and element connectivity blocks. Automatically 
%   handles point, line, surface, and volumetric elements with tailored geometric 
%   projections and boundary extraction algorithms.
%
% Input Arguments:
%   mesh - Structure representing the computational mesh, containing at least:
%     * nodes    - [N x 3] Matrix of nodal spatial coordinates (X, Y, Z).
%     * elements - Struct containing element sub-fields with connectivity matrices 
%                  and topological dimensions.

    % Validate the input mesh structure to ensure required fields exist
    validate_mesh(mesh);

    % Initialize the visualization figure with a clean white background
    figure('Name', 'Mesh Visualization', 'Color', 'w');
    hold on;
    axis equal;
    grid on;

    % Extract node coordinate matrix from the mesh data structure
    nodes = mesh.nodes;
    % Determine the spatial dimension of the domain (1D, 2D, or 3D)
    spatial_dim = get_spatial_dim(mesh, nodes);

    % Configure axes, viewing angles, and axis padding to prevent clipping
    setup_axes(spatial_dim, nodes);

    % Render elements by iterating through each element type present in the mesh
    types = fieldnames(mesh.elements);
    for t = 1:length(types)
        typeName = types{t};
        elem = mesh.elements.(typeName);
        
        % Skip element types that lack valid connectivity data
        if ~isfield(elem, 'connectivity') || isempty(elem.connectivity)
            continue;
        end
        
        % Delegate rendering to specific sub-functions based on topology
        render_element(typeName, elem, nodes, spatial_dim);
    end

    % Set a clear title showing the detected spatial dimension
    title(sprintf('Mesh Visualization (Spatial Dim: %d)', spatial_dim), ...
          'FontSize', 12, 'FontWeight', 'bold');
    hold off;
end

%% Helper Functions

function validate_mesh(mesh)
    % Ensure the input argument is a valid struct containing 'nodes' and 'elements'
    if nargin < 1 || ~isfield(mesh, 'nodes') || ~isfield(mesh, 'elements')
        error('Invalid mesh structure provided.');
    end
end

function spatial_dim = get_spatial_dim(mesh, nodes)
    % Use pre-calculated spatial_dim if available; otherwise, infer it from the columns of nodes
    if isfield(mesh, 'spatial_dim') && ~isempty(mesh.spatial_dim)
        spatial_dim = mesh.spatial_dim;
    else
        spatial_dim = size(nodes, 2);
    end
end

function setup_axes(spatial_dim, nodes)
    % Adjust the 2D/3D camera view angle depending on the spatial dimension
    if spatial_dim <= 2
        view(2); 
    else
        view(3); 
        zlabel('Z', 'FontWeight', 'bold');
    end
    
    xlabel('X', 'FontWeight', 'bold');
    ylabel('Y', 'FontWeight', 'bold');

    % Compute axis limits with 10% padding to avoid clipping boundary mesh elements
    min_coords = min(nodes, [], 1);
    max_coords = max(nodes, [], 1);
    range_coords = max_coords - min_coords;
    range_coords(range_coords == 0) = 1; % Prevent zero-division/scaling issues on flat dimensions
    
    padding = 0.1 * range_coords;
    
    xlim([min_coords(1) - padding(1), max_coords(1) + padding(1)]);
    ylim([min_coords(2) - padding(2), max_coords(2) + padding(2)]);
    
    if spatial_dim > 2 && size(nodes, 2) >= 3 && length(min_coords) >= 3
        zlim([min_coords(3) - padding(3), max_coords(3) + padding(3)]);
    end
end

function render_element(typeName, elem, nodes, spatial_dim)
    conn = elem.connectivity;
    dim = elem.dim;
    
    % Switch based on the topological dimension of the element
    switch dim
        case 1
            plot_1d(conn, nodes, typeName, spatial_dim);
        case 2
            plot_2d(conn, nodes);
        case 3
            plot_3d(conn, nodes, typeName);
        otherwise
            plot_0d(conn, nodes, spatial_dim);
    end
end

function plot_0d(conn, nodes, spatial_dim)
    % Filter out unmapped nodes (represented by zeros from read_gmsh) or out-of-bounds indices
    valid_nodes = conn(conn > 0 & conn <= size(nodes, 1));
    if ~isempty(valid_nodes)
        p = nodes(valid_nodes, :);
        % Plot 0D point markers using 2D or 3D scatter functions appropriately
        if spatial_dim <= 2 || size(nodes, 2) < 3
            scatter(p(:,1), p(:,2), 25, 'r', 'filled');
        else
            scatter3(p(:,1), p(:,2), p(:,3), 25, 'r', 'filled');
        end
    end
end

function plot_1d(conn, nodes, typeName, spatial_dim)
    % Iterate through each 1D element row
    for i = 1:size(conn, 1)
        idx = conn(i, :);
        % Filter out padding zeros and invalid node references
        idx = idx(idx > 0 & idx <= size(nodes, 1)); 
        if length(idx) >= 2
            p = nodes(idx, :);
            % Render lines matching the domain's coordinate dimensions
            if spatial_dim <= 2 || size(nodes, 2) < 3
                plot(p(:,1), p(:,2), 'b-', 'LineWidth', 1.5, ...
                     'DisplayName', sprintf('1D (%s)', typeName));
            else
                plot3(p(:,1), p(:,2), p(:,3), 'b-', 'LineWidth', 1.5, ...
                      'DisplayName', sprintf('1D (%s)', typeName));
            end
        end
    end
end

function plot_2d(conn, nodes)
    % Filter out rows containing unmapped zeros or invalid node indices to prevent patch crashes
    num_nodes = size(nodes, 1);
    valid_rows = all(conn > 0 & conn <= num_nodes, 2);
    clean_conn = conn(valid_rows, :);
    
    if isempty(clean_conn)
        return;
    end

    % Render 2D polygonal elements (e.g., triangles, quadrilaterals) using patch
    patch('Faces', clean_conn, 'Vertices', nodes, ...
          'FaceColor', [0.85 0.92 1.0], ...
          'EdgeColor', [0.15 0.15 0.4], ...
          'LineWidth', 0.8);
end

function plot_3d(conn, nodes, typeName)
    % Filter out rows containing unmapped zeros or invalid node indices
    num_nodes = size(nodes, 1);
    valid_rows = all(conn > 0 & conn <= num_nodes, 2);
    clean_conn = conn(valid_rows, :);
    
    if isempty(clean_conn)
        return;
    end

    % Render standard linear tetrahedra directly using tetramesh
    if strcmp(typeName, 'tetra4') && size(clean_conn, 2) >= 4
        tetramesh(clean_conn(:, 1:4), nodes, ...
                  'FaceColor', [0.7 0.8 1.0], ...
                  'FaceAlpha', 0.2, ...
                  'EdgeColor', [0.2 0.2 0.5], ...
                  'EdgeAlpha', 0.4, ...
                  'DisplayName', sprintf('3D (%s)', typeName));
                  
    elseif any(strcmp(typeName, {'tetra20', 'hexa8', 'prism6', 'pyramid5'}))
        if size(clean_conn, 2) >= 4
            try
                % Extract corner nodes for higher-order/volumetric elements to compute boundary faces
                tet_faces = clean_conn(:, 1:4);
                
                % Find uniquely referenced nodes and construct a local node mapping.
                % This prevents MATLAB's triangulation from issuing "Some input points are not referenced" warnings.
                used_nodes_idx = unique(tet_faces(:));
                node_map = zeros(size(nodes, 1), 1);
                node_map(used_nodes_idx) = 1:length(used_nodes_idx);
                
                compact_faces = node_map(tet_faces);
                compact_nodes = nodes(used_nodes_idx, :);
                
                % Compute boundary faces on the compact subset triangulation
                tr = triangulation(compact_faces, compact_nodes);
                boundary_faces_compact = freeBoundary(tr);
                
                % Map the compact boundary face indices back to original global node indices
                boundary_faces = used_nodes_idx(boundary_faces_compact);
                
                % Render the extracted 3D boundary mesh smoothly
                patch('Faces', boundary_faces, 'Vertices', nodes, ...
                      'FaceColor', [0.8 0.85 0.95], ...
                      'EdgeColor', [0.1 0.1 0.3], ...
                      'FaceAlpha', 0.6, ...
                      'LineWidth', 0.8, ...
                      'DisplayName', sprintf('3D Boundary (%s)', typeName));
            catch
                disp(['Note: Could not extract boundary for 3D element type ''', typeName, '''.']);
            end
        end
    else
        disp(['Note: Volumetric element type ''', typeName, ''' (3D) detected.']);
    end
end