% MATLAB equivalent of the provided Python script

%
% Open current_path in file browser
% Prompt the user to select a folder if current_path is not already defined
if ~exist('current_path', 'var') || isempty(current_path)
    current_path = uigetdir(pwd, 'Select the folder containing the image stack');
    if current_path == 0
        error('No folder selected. Operation cancelled.');
    end
end


% Function to load and process image stack from a folder


% Main Script
volume = load_images(current_path);

% Extract the mesh
threshold_value = 1;
[verts, faces] = extract_mesh(volume, threshold_value);

% Adjust scales
verts = adjust_scales(verts, current_path);

% Define plane properties
plane_position = [0, 0, 0];
plane_normal = [1, 0, 0];

% Create the figure
create_figure(verts, faces, plane_position, plane_normal);

function images = load_images(folder_path)
    % Get all .tif files in the folder
    files = dir(fullfile(folder_path, '*.tif'));

    images = [];
    for i = 1:length(files)
        file_path = fullfile(folder_path, files(i).name);
        img = imread(file_path);
        if ~isempty(img)
            images = cat(3, images, img);
        end
    end
end

% Function to extract isosurface using Marching Cubes
function [verts, faces] = extract_mesh(volume, threshold)
    % MATLAB equivalent of measure.marching_cubes
    [faces, verts] = isosurface(volume, threshold);
end

% Function to adjust scales for z-axis based on folder name
function verts = adjust_scales(verts, folder_name)
    scale_xy = 10 / 17.53;

    if contains(folder_name, 'CT')
        z_scale_factor = 4 / scale_xy;
        verts(:, 3) = verts(:, 3) * z_scale_factor;
    else
        z_scale_factor = 1 / scale_xy;
        verts(:, 3) = verts(:, 3) * z_scale_factor;
    end
end

% Function to create the figure (Plotly replacement)
function create_figure(verts, faces, plane_position, plane_normal)
    % Plot isosurface
    figure;
    p = patch('Vertices', verts, 'Faces', faces, ...
              'FaceColor', 'interp', 'EdgeColor', 'none');

    % Assign uniform color data to match the number of vertices
    numVerts = size(verts, 1);          % Number of vertices
    p.FaceVertexCData = ones(numVerts, 1); % Uniform color; can be modified as needed

    camlight;
    lighting gouraud;

    % Define slicing plane
    [xx, yy] = meshgrid( ...
        linspace(min(verts(:, 1)), max(verts(:, 1)), 50), ...
        linspace(min(verts(:, 2)), max(verts(:, 2)), 50));
    a = plane_normal(1);
    b = plane_normal(2);
    c = plane_normal(3);
    d = -dot(plane_normal, plane_position);
    zz = (-a * xx - b * yy - d) / c;

    hold on;
    surf(xx, yy, zz, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold off;

    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal;
    grid on;
end

