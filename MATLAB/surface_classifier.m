main;
function main()
    folder_path = uigetdir(pwd, 'Select the folder containing the .tif file');

    if folder_path == 0
        disp('Folder selection canceled.');
    else
        disp(['Selected folder: ', folder_path]);
    end

    output_path = "output_mesh.obj";
    if nargin < 2
        output_path = 'output_mesh.obj';
    end
    image_stack = load_image_stack(folder_path);

    mesh = create_mesh_from_stack(image_stack, folder_path);

    save_mesh(mesh, output_path);

    properties = analyze_surface_properties_from_obj(output_path);
    disp('Surface Properties from OBJ:');
    disp(properties);
end

function image_stack = load_image_stack(folder_path)
    image_files = dir(fullfile(folder_path, '*.png'));
    if isempty(image_files)
        image_files = dir(fullfile(folder_path, '*.jpg'));
    end
    if isempty(image_files)
        image_files = dir(fullfile(folder_path, '*.tif'));
    end

    image_stack = [];
    for i = 1:length(image_files)
        image = imread(fullfile(folder_path, image_files(i).name));
        if ndims(image) > 2
            image = rgb2gray(image); 
        end
        image_stack(:, :, i) = image;
    end
end

function mesh = create_mesh_from_stack(image_stack, folder_path)
    binary_stack = image_stack > 0;
    if contains(folder_path, 'CT')
        spacing = [4, 10 / 17.53, 10 / 17.53];
    else
        spacing = [1, 10 / 17.53, 10 / 17.53];
    end
    fv = isosurface(binary_stack, 0.5);
    
    scaled_vertices = fv.vertices .* spacing([3, 2, 1]); 
    mesh.vertices = scaled_vertices;
    mesh.faces = fv.faces;
end


function save_mesh(mesh, output_path)
    fid= fopen(output_path, 'w');
    fprintf(fid, 'o Mesh\n');
    for i = 1:size(mesh.vertices, 1)
        fprintf(fid, 'v %.6f %.6f %.6f\n', mesh.vertices(i, :));
    end
    for i = 1:size(mesh.faces, 1)
        fprintf(fid, 'f %d %d %d\n', mesh.faces(i, :));
    end
    fclose(fid);
end

function properties = analyze_surface_properties_from_obj(obj_path)

    [vertices, faces] = read_mesh_from_obj(obj_path);

    volume= calculate_volume(vertices, faces);
    volume_ml = abs(volume) / 1000;

    min_coords = min(vertices, [], 1);
    max_coords = max(vertices, [], 1);
    bounding_box_dims = max_coords - min_coords;

    z_dim_mm = bounding_box_dims(1);
    y_dim_mm = bounding_box_dims(2);
    x_dim_mm = bounding_box_dims(3);

    fprintf('Bounding Box Dimensions:\n');
    fprintf('  X-Dimension (mm): %.3f\n', x_dim_mm);
    fprintf('  Y-Dimension (mm): %.3f\n', y_dim_mm);
    fprintf('  Z-Dimension (mm): %.3f\n', z_dim_mm);

    edges = unique(sort([ ...
        faces(:, [1, 2]); ...
        faces(:, [2, 3]); ...
        faces(:, [3, 1])], 2), 'rows');
    edge_lengths = sqrt(sum((vertices(edges(:, 2), :) - vertices(edges(:, 1), :)).^2, 2));
    avg_edge_length = mean(edge_lengths);
    max_edge_length = max(edge_lengths);
    min_edge_length = min(edge_lengths);


    all_edges = [faces(:, [1, 2]); faces(:, [2, 3]); faces(:, [3, 1])];
   
    [~, ~, edge_indices] = unique(sort(all_edges, 2), 'rows');

    all_edges= [faces(:, [1, 2]); faces(:, [2, 3]); faces(:, [3, 1])];
    

    sorted_edges = sort(all_edges, 2);

    [edges_unique, ~, edge_indices] = unique(sorted_edges, 'rows');
    edge_face_counts = accumarray(edge_indices, 1);

    properties = struct( ...
        'surface_area_mm2', calculate_surface_area(vertices, faces), ...
        'volume_ml', volume_ml, ...
        'x_dim_mm', x_dim_mm, ...
        'y_dim_mm', y_dim_mm, ...
        'z_dim_mm', z_dim_mm, ...
        'num_triangles', size(faces, 1), ...
        'avg_edge_length_mm', avg_edge_length, ...
        'max_edge_length_mm', max_edge_length, ...
        'min_edge_length_mm', min_edge_length ...
    );
end



function [vertices, faces] = read_mesh_from_obj(file_path)
    fid = fopen(file_path, 'r');
    vertices = [];
    faces = [];
    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, 'v ')
            vertices = [vertices; sscanf(line(3:end), '%f %f %f')'];
        elseif startsWith(line, 'f ')
            faces = [faces; sscanf(line(3:end), '%d %d %d')'];
        end
    end
    fclose(fid);
end

function surface_area = calculate_surface_area(vertices, faces)
    surface_area = 0;
    for i = 1:size(faces, 1)
        v1 = vertices(faces(i, 1), :);
        v2 = vertices(faces(i, 2), :);
        v3 = vertices(faces(i, 3), :);
        surface_area = surface_area + 0.5 * norm(cross(v2 - v1, v3 - v1));
    end
end

function volume = calculate_volume(vertices, faces)
    volume = 0;
    for i = 1:size(faces, 1)
        v1 = vertices(faces(i, 1), :);
        v2 = vertices(faces(i, 2), :);
        v3 = vertices(faces(i, 3), :);
        volume = volume + dot(v1, cross(v2, v3)) / 6;
    end
end
