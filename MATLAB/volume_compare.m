clc; clear; close all;

disp('These plots are 3D so you can use the rotate tool to view them. Please upload the .obj file from surfaceclassifier')
[filename, pathname] = uigetfile({'*.obj'}, 'Select an OBJ file');
if isequal(filename, 0)
    disp('User canceled file selection.');
    return;
end
filepath = fullfile(pathname, filename);


[obj.vertices, obj.faces] = readObj(filepath);

obj.vertices = obj.vertices * 0.1; 

imported_volume = computeMeshVolume(obj.vertices, obj.faces);

minCoords = min(obj.vertices);
maxCoords = max(obj.vertices);
objSize = maxCoords - minCoords; 
objRadius = max(objSize) / 2;


volumes = struct( ...
    'tennis_ball', 4/3 * pi * (6.7/2)^3, ...      % Tennis ball (6.7 cm diameter)
    'football', 4/3 * pi * (11/2)^2 * (28/2), ... % Football (11 cm minor axis, 28 cm major axis)
    'basketball', 4/3 * pi * (24.6/2)^3, ...      % Basketball (24.6 cm diameter)
    'abdomen', 4/3 * pi * (20/2) * (35/2) * (25/2) ... % Human abdomen (approximate ellipsoid)
);

fprintf('Compared to:\n');
fprintf('  - Tennis Ball: %.2fx\n', imported_volume / volumes.tennis_ball);
fprintf('  - Football: %.2fx\n', imported_volume / volumes.football);
fprintf('  - Basketball: %.2fx\n', imported_volume / volumes.basketball);
fprintf('  - Human Abdomen: %.2fx\n', imported_volume / volumes.abdomen);


% Tennis Ball
figure;
hold on;
tennisRadius = (3 * volumes.tennis_ball / (4 * pi))^(1/3);
plotSphere(volumes.tennis_ball, [0, 0, 0], 'g');
plotImportedObject(obj, [0, 0, 0], tennisRadius);
title('Imported Object vs. Tennis Ball');
hold off;

% Basketball
figure;
hold on;
basketballRadius = (3 * volumes.basketball / (4 * pi))^(1/3);
plotSphere(volumes.basketball, [0, 0, 0], 'b');
plotImportedObject(obj, [0, 0, 0], basketballRadius);
title('Imported Object vs. Basketball');
hold off;

% Football
figure;
hold on;
footballRadius = (volumes.football / (4/3 * pi))^(1/3); 
plotEllipsoid(volumes.football, [0, 0, 0], 'm');
plotImportedObject(obj, [0, 0, 0], footballRadius);
title('Imported Object vs. Football');
hold off;


% Human Full-Body Outline with Imported Object in Abdomen
figure;
hold on;
plotFullHumanOutline();
plotImportedObjectInAbdomen(obj);
title('Imported Object in Human Body');
hold off;


function [vertices, faces] = readObj(filename)
    fid = fopen(filename);
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

function volume = computeMeshVolume(vertices, faces)
    volume = 0;
    for i = 1:size(faces, 1)
        v1 = vertices(faces(i, 1), :);
        v2 = vertices(faces(i, 2), :);
        v3 = vertices(faces(i, 3), :);
        volume = volume + dot(v1, cross(v2, v3)) / 6;
    end
    volume = abs(volume);
end

function plotSphere(volume, center, color)
    radius = (3 * volume / (4 * pi))^(1/3);
    [X, Y, Z] = sphere;
    X = X * radius + center(1);
    Y = Y * radius + center(2);
    Z = Z * radius + center(3);
    surf(X, Y, Z, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    axis equal;
    xlabel('X (cm)'); ylabel('Y (cm)'); zlabel('Z (cm)');
    camlight; lighting gouraud;
end

function plotEllipsoid(volume, center, color)
    [X, Y, Z] = ellipsoid(center(1), center(2), center(3), 1, 1.5, 1.2, 30);
    scaleFactor = (volume / (4/3 * pi))^(1/3);
    X = X * scaleFactor;
    Y = Y * scaleFactor;
    Z = Z * scaleFactor;
    surf(X, Y, Z, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    axis equal;
    xlabel('X (cm)'); ylabel('Y (cm)'); zlabel('Z (cm)');
    camlight; lighting gouraud;
end

function plotImportedObject(obj, referenceCenter, referenceRadius)
    minCoords = min(obj.vertices);
    maxCoords = max(obj.vertices);
    objSize = maxCoords - minCoords;
    
    objCenter = (minCoords + maxCoords) / 2;

    objRadius = max(objSize) / 2;

    offset = referenceRadius + objRadius;
    newPosition = referenceCenter + [offset - objCenter(1), -objCenter(2), -objCenter(3)];

    translated_vertices = obj.vertices + newPosition;
    patch('Faces', obj.faces, 'Vertices', translated_vertices, ...
          'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    axis equal;
end



function plotFullHumanOutline()
    % Define proportions (height, widths, radii)
    human_height = 175;
    shoulder_width = 0.5 * human_height;
    waist_width = 0.3 * human_height;
    leg_spacing = 0.2 * human_height;
    leg_length = 0.45 * human_height;
    torso_length = 0.5 * human_height;
    head_radius = 0.15 * human_height;

    head_center = [0, human_height - head_radius];
    shoulder_y = human_height - 2 * head_radius;
    waist_y = shoulder_y - torso_length;
    leg_start_y = waist_y - 0.05 * human_height;
    left_leg_x = [-waist_width/2, -waist_width/2];
    right_leg_x = [waist_width/2, waist_width/2];
    leg_y = [leg_start_y, leg_start_y - leg_length];

    arm_length = 0.4 * human_height;
    left_arm_x = [-shoulder_width/2, -shoulder_width/2 - 15];
    right_arm_x = [shoulder_width/2, shoulder_width/2 + 15];
    arm_y = [shoulder_y, shoulder_y - arm_length];

    theta = linspace(0, 2*pi, 100);
    plot(head_center(1) + head_radius*cos(theta), head_center(2) + head_radius*sin(theta), 'k', 'LineWidth', 4);
    hold on;
    plot([-shoulder_width/2, -waist_width/2], [shoulder_y, waist_y], 'k', 'LineWidth', 4);
    plot([shoulder_width/2, waist_width/2], [shoulder_y, waist_y], 'k', 'LineWidth', 4);
    plot(left_leg_x, leg_y, 'k', 'LineWidth', 4);
    plot(right_leg_x, leg_y, 'k', 'LineWidth', 4);
    plot(left_arm_x, arm_y, 'k', 'LineWidth', 4);
    plot(right_arm_x, arm_y, 'k', 'LineWidth', 4);
    axis equal;
    xlabel('X (cm)'); ylabel('Y (cm)');
    title('2D Human Outline');
    hold off;
end




function plotImportedObjectInAbdomen(obj)
    theta = pi/2;
    R = [cos(theta)  0  sin(theta);
         0           1  0;
        -sin(theta)  0  cos(theta)];
    objCenter = mean(obj.vertices, 1);
    
    obj.vertices = (obj.vertices - objCenter) * R' + objCenter;

    minCoords = min(obj.vertices);
    maxCoords = max(obj.vertices);
    objSize = maxCoords - minCoords;
    
    objCenter = mean(obj.vertices, 1);

    abdomenCenter = [0, 80, 0]; 
    abdomenRadius = 15; 

    offset = abdomenRadius;
    newPosition = abdomenCenter - objCenter + [0, -offset, 0];

    translated_vertices = obj.vertices + newPosition;

    patch('Faces', obj.faces, 'Vertices', translated_vertices, ...
          'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);


    view(3);
    axis equal;
    xlabel('X (cm)');
    ylabel('Y (cm)');
    zlabel('Z (cm)');
    title('Imported 3D Object in 2D Human Outline');

    warning_text = 'Note: The location shown is for visualization only.';
    text(0, -20, warning_text, 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

