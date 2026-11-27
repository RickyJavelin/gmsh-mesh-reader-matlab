function mesh = read_gmsh(filename)
% READ_GMSH Parse Gmsh mesh files into a structured format for finite element analysis.
%
% Syntax:
%   mesh = read_gmsh(filename)
%
% Description:
%   Reads an ASCII Gmsh mesh file (supporting versions 2.x and 4.x), extracts 
%   nodal coordinates, global node identifiers, element connectivity matrices, 
%   and determines the topological and geometric spatial dimensions.
%
% Input Arguments:
%   filename - Character vector or string specifying the relative or absolute 
%              path to the Gmsh (.msh) file.
%
% Output Arguments:
%   mesh - Structure containing the following fields:
%     * node_tags   - [N x 1] Column vector of global node identifiers.
%     * nodes       - [N x 3] Matrix of nodal coordinates (X, Y, Z); 
%                     the third coordinate is set to zero for two-dimensional domains.
%     * elements    - Struct containing fields corresponding to each element type 
%                     present in the mesh. Each element substructure comprises:
%                     - tags         - [Ne x 1] Vector of element identifiers.
%                     - dim          - Topological dimension of the element.
%                     - connectivity - [Ne x Nv] Matrix of nodal indices defining 
%                                      the element connectivity mapped to local indices.
%     * spatial_dim - Geometric spatial dimension of the computational domain.

    % Load the entire file into memory as a string array for efficient vectorised parsing
    lines = readlines(filename);
    
    % Identify the Gmsh file version from the $MeshFormat section header
    format_start = find(lines == "$MeshFormat", 1);
    if ~isempty(format_start)
        format_header = sscanf(lines(format_start + 1), '%f');
        file_version = format_header(1);
    else
        warning('No mesh format specification found; assuming default compatibility.');
    end
    
    % Locate the beginning of the node coordinate data block
    node_start = find(lines == "$Nodes", 1);
    if isempty(node_start)
        error('Invalid Gmsh file structure: $Nodes section not found.');
    end
    
    % Locate the beginning of the element connectivity data block
    elem_start = find(lines == "$Elements", 1);
    if isempty(elem_start)
        error('Invalid Gmsh file structure: $Elements section not found.');
    end
    
    %% Parse Node Data
    if floor(file_version) == 2
        % Version 2.x: Sequential block format
        node_header = sscanf(lines(node_start + 1), '%f');
        numNodesTotal = node_header(1);
        
        node_data_text = lines(node_start + 2 : node_start + 1 + numNodesTotal);
        node_data = sscanf(strjoin(node_data_text, '\n'), '%f');
        node_matrix = reshape(node_data, 4, [])'; 
        
        node_tags = node_matrix(:, 1);
        node_coords = node_matrix(:, 2:4);
    else
        % Version 4.x+: Entity-based block structured format
        node_header = sscanf(lines(node_start + 1), '%f');
        numEntityBlocks = node_header(1);
        numNodesTotal   = node_header(2);
        
        node_tags = zeros(numNodesTotal, 1);
        node_coords = zeros(numNodesTotal, 3);
        
        curr_line = node_start + 2;
        node_idx = 1;
        
        for b = 1:numEntityBlocks
            block_header = sscanf(lines(curr_line), '%f');
            numNodesInBlock = block_header(4);
            curr_line = curr_line + 1;
            
            tags_block = sscanf(strjoin(lines(curr_line : curr_line + numNodesInBlock - 1), '\n'), '%f');
            curr_line = curr_line + numNodesInBlock;
            
            coords_text = lines(curr_line : curr_line + numNodesInBlock - 1);
            coords_block = reshape(sscanf(strjoin(coords_text, '\n'), '%f'), 3, [])'; 
            curr_line = curr_line + numNodesInBlock;
            
            range = node_idx : node_idx + numNodesInBlock - 1;
            node_tags(range) = tags_block;
            node_coords(range, :) = coords_block;
            node_idx = node_idx + numNodesInBlock;
        end
    end
    
    % Sort nodes by their global identifiers to ensure monotonic indexing consistency
    [mesh.node_tags, sort_idx] = sort(node_tags);
    mesh.nodes = node_coords(sort_idx, :);
    
    %% Parse Element Data
    elements = struct();
    if floor(file_version) == 2
        % Version 2.x element parsing workflow
        elem_header = sscanf(lines(elem_start + 1), '%f');
        numElementsTotal = elem_header(1);
        
        elem_text = lines(elem_start + 2 : elem_start + 1 + numElementsTotal);
        
        % Extract metadata columns (element number, element type, number of tags)
        header_data = cell2mat(cellfun(@(x) sscanf(x, '%f', 3)', elem_text, 'UniformOutput', false));

        element_types = header_data(:, 2);
        num_tags_col = header_data(:, 3);
        unique_types = unique(element_types)';
        
        for t = 1:length(unique_types)
            elementType = unique_types(t);
            mask = (element_types == elementType);
            
            sub_text = elem_text(mask);
            sub_num_tags = num_tags_col(mask);
            num_tags_val = sub_num_tags(1);
            
            sub_block_str = strjoin(sub_text, '\n');
            sub_vals = sscanf(sub_block_str, '%f');
            
            items_per_line = length(sub_vals) / length(sub_text);
            data_matrix = reshape(sub_vals, items_per_line, [])';
            
            elem_tags_type = data_matrix(:, 1);
            node_start_idx = 4 + num_tags_val;
            raw_nodes = data_matrix(:, node_start_idx:end);
            
            [typeName, elemDim] = getElementTypeAndDim(elementType);
            elements.(typeName).tags = elem_tags_type;
            elements.(typeName).dim = elemDim;
            
            % Map global node tags in element connectivity to sorted local node indices
            [~, connectivity] = ismember(raw_nodes, mesh.node_tags);
            elements.(typeName).connectivity = connectivity;
        end
        
    else
        % Version 4.x+ element parsing via entity blocks workflow
        elem_header = sscanf(lines(elem_start + 1), '%f');
        numElemEntityBlocks = elem_header(1);
        
        curr_line = elem_start + 2;
        temp_data = struct();
        
        for b = 1:numElemEntityBlocks
            block_header = sscanf(lines(curr_line), '%f');
            elementType        = block_header(3);
            numElementsInBlock = block_header(4);
            curr_line = curr_line + 1;
            
            elem_lines = lines(curr_line : curr_line + numElementsInBlock - 1);
            curr_line = curr_line + numElementsInBlock;
            
            block_data = sscanf(strjoin(elem_lines, '\n'), '%f');
            
            col_len = length(block_data) / numElementsInBlock;
            data_matrix = reshape(block_data, col_len, numElementsInBlock)';
            
            elem_tags = data_matrix(:, 1);
            raw_nodes = data_matrix(:, 3:end);
            
            % Map global node tags to local indices
            [~, connectivity] = ismember(raw_nodes, mesh.node_tags);
            
            [typeName, elemDim] = getElementTypeAndDim(elementType);
            if ~isfield(temp_data, typeName)
                temp_data.(typeName).tags = {};
                temp_data.(typeName).connectivity = {};
                temp_data.(typeName).dim = elemDim;
            end
            temp_data.(typeName).tags{end+1, 1} = elem_tags;
            temp_data.(typeName).connectivity{end+1, 1} = connectivity;
        end
        
        types = fieldnames(temp_data);
        for t = 1:length(types)
            typeName = types{t};
            elements.(typeName).tags = cell2mat(temp_data.(typeName).tags);
            elements.(typeName).connectivity = cell2mat(temp_data.(typeName).connectivity);
            elements.(typeName).dim = temp_data.(typeName).dim;
        end
    end
    
    mesh.elements = elements;

    %% Determine Effective Spatial Dimension
    max_coord_col = size(mesh.nodes, 2);
    if max_coord_col == 1
        mesh.spatial_dim = 1;
    elseif max_coord_col == 2
        if all(mesh.nodes(:, 2) == 0)
            mesh.spatial_dim = 1;
        else
            mesh.spatial_dim = 2;
        end
    else
        if all(mesh.nodes(:, 3) == 0)
            if all(mesh.nodes(:, 2) == 0)
                mesh.spatial_dim = 1;
            else
                mesh.spatial_dim = 2;
            end
        else
            mesh.spatial_dim = 3;
        end
    end

end

function [name, dim] = getElementTypeAndDim(typeNum)
    % GETELEMENTTYPEANDDIM Map Gmsh numerical element type identifiers to standard nomenclature and topological dimensions.
    switch typeNum
        case 1, name = 'line2';     dim = 1;
        case 2, name = 'triangle3'; dim = 2;
        case 3, name = 'quad4';     dim = 2;
        case 4, name = 'tetra4';    dim = 3;
        case 5, name = 'hexa8';     dim = 3;
        case 6, name = 'prism6';    dim = 3;
        case 7, name = 'pyramid5';  dim = 3;
        case 8, name = 'line3';     dim = 1;
        case 9, name = 'triangle6'; dim = 2;
        case 15, name = 'node1';    dim = 0;
        case 29, name = 'tetra20';  dim = 3;
        otherwise
            name = sprintf('type_%d', typeNum);
            dim = 2;
    end
end