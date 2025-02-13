function texture_analysis_main()
% Main function to process CT or MRI images from a folder, compute texture
% features, save results to CSV, and display the images with analysis.
%
% Change the image_directory and modality ('CT' or 'MRI') as needed.

    %% Set Parameters
    image_directory = '/Users/colehanan/Desktop/Important_AMOS/slicedUpImage500/group_4/CT_2';
    modality = 'CT';  % Change to 'MRI' if processing MRI images

    % For CT images, create (or specify) a directory to save converted images.
    if strcmpi(modality, 'CT')
        converted_directory = fullfile(image_directory, 'CT_2_converted');
        if ~exist(converted_directory, 'dir')
            mkdir(converted_directory);
        end
        output_csv = fullfile(converted_directory, 'texture_analysis_results.csv');
    else
        % For MRI, results will be saved in the same directory.
        converted_directory = [];
        output_csv = fullfile(image_directory, 'texture_analysis_results.csv');
    end

    %% Analyze the image stack
    all_features = analyze_image_stack(image_directory, modality, output_csv, converted_directory);

    %% Display results
    if strcmpi(modality, 'CT')
        show_ct_images_with_analysis(converted_directory, all_features);
    elseif strcmpi(modality, 'MRI')
        show_mri_images_with_analysis(image_directory, all_features);
    end
end

%% ====================================================
function output_path = ensure_ct_image(input_path, output_path)
% Reads a CT image and ensures it is in one of the allowed types.
% If not (e.g. if it is 64-bit or another unsupported type), it is cast
% to single (32-bit float). The result is then saved.
%
% Parameters:
%   input_path:  full path to the original image
%   output_path: full path to save the converted image
%
% Returns:
%   output_path (nonempty if successful, empty if an error occurred)

    try
        image = imread(input_path);
    catch ME
        fprintf('Error reading CT image %s: %s\n', input_path, ME.message);
        output_path = '';
        return;
    end

    % Allowed classes for TIFF writing: uint8, int8, uint16, int16, single
    if ~(isa(image, 'uint8') || isa(image, 'int8') || ...
         isa(image, 'uint16') || isa(image, 'int16') || isa(image, 'single'))
        % Convert to single (32-bit float)
        image = im2single(image);
        fprintf('Casting CT image %s to single (32-bit float).\n', input_path);
    end

    % Write image to output_path as TIFF
    try
        imwrite(image, output_path);
        fprintf('Saved CT converted image to %s\n', output_path);
    catch ME
        fprintf('Error writing CT image to %s: %s\n', output_path, ME.message);
        output_path = '';
        return;
    end
end

%% ====================================================
function features = analyze_ct_texture(image_path)
% Analyzes the texture of a CT image.
% The image is converted to 8-bit (if necessary) and then used to compute
% the gray-level co-occurrence matrix (GLCM) and texture properties.
%
% Parameters:
%   image_path: full path to the CT image
%
% Returns:
%   features: a structure with texture properties (contrast, correlation,
%             energy, homogeneity, mean, std, median)

    try
        image = imread(image_path);
    catch ME
        fprintf('Error reading CT image %s: %s\n', image_path, ME.message);
        features = [];
        return;
    end

    if max(image(:)) == 0
        fprintf('Skipping CT image %s because it is a black frame.\n', image_path);
        features = [];
        return;
    end

    % Convert image to uint8 for graycomatrix (if not already)
    image_uint8 = im2uint8(image);

    % Define offsets for a distance of 5 pixels at angles: 0°, 45°, 90°, 135°
    offsets = [0 5; -5 5; -5 0; -5 -5];
    glcm = graycomatrix(image_uint8, 'Offset', offsets, 'Symmetric', true, 'NumLevels', 256);
    stats = graycoprops(glcm, {'Contrast','Correlation','Energy','Homogeneity'});

    % Average across all offsets
    features.contrast = mean(stats.Contrast);
    features.correlation = mean(stats.Correlation);
    features.energy = mean(stats.Energy);
    features.homogeneity = mean(stats.Homogeneity);

    features.mean = mean(image_uint8(:));
    features.std = std(double(image_uint8(:)));
    features.median = median(image_uint8(:));
end

%% ====================================================
function features = analyze_mri_texture(image_path)
% Analyzes the texture of an MRI image.
% The image is read (even if 64-bit), then min–max normalized to 0–255 and
% converted to 8-bit before computing the GLCM and texture properties.
%
% Parameters:
%   image_path: full path to the MRI image
%
% Returns:
%   features: a structure with texture properties

    try
        image = imread(image_path);
    catch ME
        fprintf('Error reading MRI image %s: %s\n', image_path, ME.message);
        features = [];
        return;
    end

    if max(image(:)) == 0
        fprintf('Skipping MRI image %s because it is a black frame.\n', image_path);
        features = [];
        return;
    end

    image = double(image);
    % Min–max normalization to [0,1]
    img_min = min(image(:));
    img_max = max(image(:));
    if (img_max - img_min) ~= 0
        norm_image = (image - img_min) / (img_max - img_min);
    else
        norm_image = image - img_min;
    end
    image_uint8 = uint8(norm_image * 255);

    offsets = [0 5; -5 5; -5 0; -5 -5];
    glcm = graycomatrix(image_uint8, 'Offset', offsets, 'Symmetric', true, 'NumLevels', 256);
    stats = graycoprops(glcm, {'Contrast','Correlation','Energy','Homogeneity'});

    features.contrast = mean(stats.Contrast);
    features.correlation = mean(stats.Correlation);
    features.energy = mean(stats.Energy);
    features.homogeneity = mean(stats.Homogeneity);

    features.mean = mean(image_uint8(:));
    features.std = std(double(image_uint8(:)));
    features.median = median(image_uint8(:));
end

%% ====================================================
function all_features = analyze_image_stack(image_directory, modality, output_csv, converted_dir)
% Analyzes a stack of images in a given directory and saves the results
% to a CSV file. For CT images, each file is first “converted” if needed.
%
% Parameters:
%   image_directory: folder containing the TIFF images
%   modality: 'CT' or 'MRI'
%   output_csv: path to save the CSV file
%   converted_dir: (for CT) folder to save converted images
%
% Returns:
%   all_features: a structure array with one element per processed image

    % Get list of .tif files
    files = dir(fullfile(image_directory, '*.tif'));
    filenames = {files.name};

    % Sort filenames based on a numeric token (e.g., '_0001.tif')
    numbers = zeros(length(filenames),1);
    for i = 1:length(filenames)
        token = regexp(filenames{i}, '_([\d]{4})\.tif$', 'tokens');
        if ~isempty(token)
            numbers(i) = str2double(token{1}{1});
        else
            numbers(i) = 0;
        end
    end
    [~, idx] = sort(numbers);
    filenames = filenames(idx);

    % Initialize structure array for features
    all_features = struct('filename', {}, 'contrast', {}, 'correlation', {}, ...
        'energy', {}, 'homogeneity', {}, 'mean', {}, 'std', {}, 'median', {});

    if strcmpi(modality, 'CT')
        if nargin < 4 || isempty(converted_dir)
            converted_dir = fullfile(image_directory, 'CT_converted');
            if ~exist(converted_dir, 'dir')
                mkdir(converted_dir);
            end
        end
        for i = 1:length(filenames)
            original_path = fullfile(image_directory, filenames{i});
            converted_path = fullfile(converted_dir, filenames{i});
            ct_path = ensure_ct_image(original_path, converted_path);
            if isempty(ct_path)
                fprintf('Skipping %s due to conversion error.\n', original_path);
                continue;
            end
            features = analyze_ct_texture(ct_path);
            if ~isempty(features)
                newFeature.filename = filenames{i};
                newFeature.contrast = features.contrast;
                newFeature.correlation = features.correlation;
                newFeature.energy = features.energy;
                newFeature.homogeneity = features.homogeneity;
                newFeature.mean = features.mean;
                newFeature.std = features.std;
                newFeature.median = features.median;
                all_features(end+1) = newFeature; %#ok<AGROW>
            end
        end
    elseif strcmpi(modality, 'MRI')
        for i = 1:length(filenames)
            image_path = fullfile(image_directory, filenames{i});
            features = analyze_mri_texture(image_path);
            if ~isempty(features)
                newFeature.filename = filenames{i};
                newFeature.contrast = features.contrast;
                newFeature.correlation = features.correlation;
                newFeature.energy = features.energy;
                newFeature.homogeneity = features.homogeneity;
                newFeature.mean = features.mean;
                newFeature.std = features.std;
                newFeature.median = features.median;
                all_features(end+1) = newFeature; %#ok<AGROW>
            end
        end
    else
        disp('Unsupported modality. Use "CT" or "MRI".');
        all_features = [];
        return;
    end

    % Convert the structure array to a table and write CSV.
    T = struct2table(all_features);
    writetable(T, output_csv);
    fprintf('Texture analysis results saved to %s\n', output_csv);
end

%% ====================================================
function show_ct_images_with_analysis(converted_dir, all_features)
% Displays CT images from the converted folder with their texture analysis.
%
% Parameters:
%   converted_dir: directory containing the converted CT images
%   all_features: structure array of texture features

    for i = 1:length(all_features)
        image_file = all_features(i).filename;
        image_path = fullfile(converted_dir, image_file);
        try
            image = imread(image_path);
        catch ME
            fprintf('Error reading CT image %s for display: %s\n', image_path, ME.message);
            continue;
        end

        figure;
        imshow(image, []);
        title(sprintf('%s\nContrast: %.2f, Correlation: %.2f,\nEnergy: %.2f, Homogeneity: %.2f',...
            image_file, all_features(i).contrast, all_features(i).correlation, ...
            all_features(i).energy, all_features(i).homogeneity));
    end
end

%% ====================================================
function show_mri_images_with_analysis(image_directory, all_features)
% Displays MRI images with their texture analysis.
%
% Parameters:
%   image_directory: folder containing the MRI images
%   all_features: structure array of texture features

    for i = 1:length(all_features)
        image_file = all_features(i).filename;
        image_path = fullfile(image_directory, image_file);
        try
            image = imread(image_path);
            image = double(image);
            img_min = min(image(:));
            img_max = max(image(:));
            if (img_max - img_min) ~= 0
                norm_image = (image - img_min) / (img_max - img_min);
            else
                norm_image = image - img_min;
            end
            image_uint8 = uint8(norm_image * 255);
        catch ME
            fprintf('Error reading MRI image %s for display: %s\n', image_path, ME.message);
            continue;
        end

        figure;
        imshow(image_uint8, []);
        title(sprintf('%s\nContrast: %.2f, Correlation: %.2f,\nEnergy: %.2f, Homogeneity: %.2f',...
            image_file, all_features(i).contrast, all_features(i).correlation, ...
            all_features(i).energy, all_features(i).homogeneity));
    end
end
