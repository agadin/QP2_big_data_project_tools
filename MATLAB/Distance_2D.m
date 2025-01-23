function Distance_2D
    % In the file window that first appears select the folder in which the
    % images are stored in.
    % In the window with the images, hit the Measure Distance button to
    % take a measurement.
    folder_path = uigetdir('', 'Select a Folder Containing .tif Files');
    if folder_path == 0
        disp('No folder selected. Exiting.');
        return;
    end

    tif_files = dir(fullfile(folder_path, '*.tif'));
    if isempty(tif_files)
        disp(['No .tif files found in ', folder_path]);
        return;
    end

    current_index = 1;
    fig = figure('Name', 'Measure Distance on Image with File Navigation', ...
                 'NumberTitle', 'off', 'Position', [100, 100, 800, 800]);

    ax = axes('Parent', fig, 'Position', [0.1, 0.3, 0.8, 0.6]);
    img = imread(fullfile(folder_path, tif_files(current_index).name));
    imshow(img, [], 'Parent', ax);
    title(ax, 'Click two points to measure distance');

    slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', length(tif_files), ...
                       'Value', current_index, 'SliderStep', [1/(length(tif_files)-1), 1/(length(tif_files)-1)], ...
                       'Position', [100, 50, 600, 20], 'Parent', fig);
    addlistener(slider, 'Value', 'PostSet', @(src, event) update_image());

    uicontrol('Style', 'pushbutton', 'String', 'Measure Distance', ...
              'Position', [350, 10, 100, 30], 'Callback', @measure_distance);

    function update_image
        current_index = round(slider.Value);
        img = imread(fullfile(folder_path, tif_files(current_index).name));
        imshow(img, [], 'Parent', ax);
        title(ax, 'Click two points to measure distance');
    end
    function measure_distance(~, ~)
        disp('Click two points on the image...');
        [x, y] = ginput(2); 
        hold on;
        plot(x, y, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
        line(x, y, 'Color', 'r', 'LineWidth', 2);
        hold off;

        pixel_distance = sqrt((x(2) - x(1))^2 + (y(2) - y(1))^2);
        mm_distance = pixel_distance * (10 / 17.52); 

        disp(['Distance: ', num2str(pixel_distance, '%.2f'), ' pixels (', ...
              num2str(mm_distance, '%.2f'), ' mm)']);
    end
end