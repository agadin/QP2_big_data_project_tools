function all_features = textureAnalysis(image_directory)
% textureAnalysis performs texture analysis on PNG CT images.
%
% USAGE:
%   all_features = textureAnalysis(image_directory)
%
% INPUT:
%   image_directory - (string) Folder containing your PNG images.
%                     If not provided, a folder selection dialog appears.
%
% OUTPUT:
%   all_features - A structure containing texture features for each image.
%
% The function:
%   - Creates a subfolder 'CT_converted' for converted images.
%   - Skips images that are completely black.
%   - Computes texture features (Contrast, Correlation, Energy, Homogeneity)
%     from the Gray-Level Co-occurrence Matrix (GLCM) and basic statistics
%     (Mean, Std, Median).
%   - Saves the results to a CSV file in the converted folder.
%   - Displays each image with an overlayed title showing key metrics.

    if nargin < 1 || isempty(image_directory)
        image_directory = uigetdir(pwd, 'Select folder containing PNG images');
        if isequal(image_directory, 0)
            error('No directory selected.');
        end
    end

    % Create a subfolder for converted images.
    converted_dir = fullfile(image_directory, 'CT_converted');
    if ~exist(converted_dir, 'dir')
        mkdir(converted_dir);
    end

    % Define the CSV output file.
    output_csv = fullfile(converted_dir, 'texture_analysis_results.csv');

    % Analyze the image stack.
    all_features = analyze_image_stack(image_directory, output_csv, converted_dir);

    % Display images with the computed texture analysis.
    show_ct_images_with_analysis(converted_dir, all_features);
end

%% Subfunction: Ensure CT image is in a supported format
function output_path = ensure_ct_image(input_path, output_path)
% Reads a PNG image and ensures its pixel type is allowed.
% If not (i.e. not uint8, int8, uint16, int16, or single),
% the image is converted to single precision.
% The resulting image is saved to output_path.

    try
        image = imread(input_path);
    catch ME
        fprintf('Error reading image %s: %s\n', input_path, ME.message);
        output_path = '';
        return;
    end

    % Allowed types: uint8, int8, uint16, int16, single.
    allowed_types = {'uint8','int8','uint16','int16','single'};
    if ~any(strcmp(class(image), allowed_types))
        image = im2single(image);
        fprintf('Casting image %s to single precision.\n', input_path);
    end

    try
        imwrite(image, output_path);
        fprintf('Saved converted image to %s\n', output_path);
    catch ME
        fprintf('Error writing image to %s: %s\n', output_path, ME.message);
        output_path = '';
        return;
    end
end

%% Subfunction: Analyze CT texture from one image
function features = analyze_ct_texture(image_path)
% Reads an image, converts it to 8-bit, computes the GLCM, and extracts texture
% features along with basic statistics. Completely black images are skipped.

    try
        image = imread(image_path);
    catch ME
        fprintf('Error reading image %s: %s\n', image_path, ME.message);
        features = [];
        return;
    end

    % If image is color, convert to grayscale.
    if size(image,3) == 3
        image = rgb2gray(image);
    end

    % Check if the image is completely black.
    if max(image(:)) == 0
        fprintf('Skipping image %s because it is a black frame.\n', image_path);
        features = [];
        return;
    end

    % Convert image to uint8 if not already.
    if ~isa(image, 'uint8')
        image = im2uint8(image);
    end

    % Define offsets corresponding to distance 5 and angles: 0, 45, 90, 135 degrees.
    % Offsets are defined as [row, col] differences.
    offset1 = [0, 5];                     % Angle 0°
    offset2 = [-round(5*sin(pi/4)), round(5*cos(pi/4))]; % Angle 45° (approx. [-4, 4])
    offset3 = [-5, 0];                    % Angle 90°
    offset4 = [-round(5*sin(3*pi/4)), -round(5*cos(3*pi/4))]; % Angle 135° (approx. [-4, -4])
    offsets = [offset1; offset2; offset3; offset4];

    % Compute the Gray-Level Co-occurrence Matrix (GLCM).
    glcm = graycomatrix(image, 'Offset', offsets, 'Symmetric', true, ...
                         'NumLevels', 256, 'GrayLimits', [0 255]);
    % Compute texture properties.
    stats = graycoprops(glcm, {'Contrast','Correlation','Energy','Homogeneity'});
    % Use properties from the first offset (to mimic the Python [0,0] selection).
    contrast    = stats.Contrast(1);
    correlation = stats.Correlation(1);
    energy      = stats.Energy(1);
    homogeneity = stats.Homogeneity(1);

    % Compute basic statistics.
    mean_val   = mean(image(:));
    std_val    = std(double(image(:)));
    median_val = median(image(:));

    % Return the features in a structure.
    features = struct('contrast', contrast, ...
                      'correlation', correlation, ...
                      'energy', energy, ...
                      'homogeneity', homogeneity, ...
                      'mean', mean_val, ...
                      'std', std_val, ...
                      'median', median_val);
end

%% Subfunction: Analyze a stack of images in a directory
function all_features = analyze_image_stack(image_directory, output_csv, converted_dir)
% Processes all PNG images in the specified directory.
% Each image is converted (if necessary) and analyzed.
% Results are saved to a CSV file.

    all_features = struct();

    % List all PNG files in the directory.
    files = dir(fullfile(image_directory, '*.png'));

    % Sort files based on a numerical token in the filename.
    % The regular expression looks for an underscore followed by 4 digits before ".png".
    numbers = zeros(length(files),1);
    for i = 1:length(files)
        token = regexp(files(i).name, '_(\d{4})\.png$', 'tokens');
        if ~isempty(token)
            numbers(i) = str2double(token{1}{1});
        else
            numbers(i) = 0;
        end
    end
    [~, sortIdx] = sort(numbers);
    files = files(sortIdx);

    % Ensure the converted images directory exists.
    if ~exist(converted_dir, 'dir')
        mkdir(converted_dir);
    end

    % Process each image.
    for i = 1:length(files)
        image_file = files(i).name;
        original_path = fullfile(image_directory, image_file);
        converted_path = fullfile(converted_dir, image_file);

        ct_path = ensure_ct_image(original_path, converted_path);
        if isempty(ct_path)
            fprintf('Skipping %s due to conversion error.\n', original_path);
            continue;
        end

        features = analyze_ct_texture(ct_path);
        if ~isempty(features)
            % Use a valid fieldname (replace any invalid characters).
            fieldName = matlab.lang.makeValidName(image_file);
            all_features.(fieldName) = features;
        end
    end

    % Convert the results into a table and write to a CSV file.
    fileNames = fieldnames(all_features);
    if isempty(fileNames)
        fprintf('No features to save.\n');
        return;
    end

    % Preallocate arrays for table columns.
    contrastArr = [];
    correlationArr = [];
    energyArr = [];
    homogeneityArr = [];
    meanArr = [];
    stdArr = [];
    medianArr = [];
    names = cell(length(fileNames), 1);

    for i = 1:length(fileNames)
        fname = fileNames{i};
        feat = all_features.(fname);
        contrastArr(i,1)    = feat.contrast;
        correlationArr(i,1) = feat.correlation;
        energyArr(i,1)      = feat.energy;
        homogeneityArr(i,1) = feat.homogeneity;
        meanArr(i,1)        = feat.mean;
        stdArr(i,1)         = feat.std;
        medianArr(i,1)      = feat.median;
        names{i}            = fname;
    end

    T = table(names, contrastArr, correlationArr, energyArr, homogeneityArr, ...
              meanArr, stdArr, medianArr, ...
              'VariableNames',{'Filename','Contrast','Correlation','Energy','Homogeneity','Mean','Std','Median'});
    writetable(T, output_csv);
    fprintf('Texture analysis results saved to %s\n', output_csv);
end

%% Subfunction: Display CT images with analysis
function show_ct_images_with_analysis(converted_dir, all_features)
% Displays each converted image along with its computed texture features.
    fileNames = fieldnames(all_features);
    for i = 1:length(fileNames)
        image_file = fileNames{i};
        image_path = fullfile(converted_dir, image_file);
        try
            image = imread(image_path);
        catch ME
            fprintf('Error reading image %s for display: %s\n', image_path, ME.message);
            continue;
        end

        feat = all_features.(image_file);
        figure;
        imshow(image);
        title(sprintf('%s\nContrast: %.2f, Correlation: %.2f, Energy: %.2f, Homogeneity: %.2f', ...
            image_file, feat.contrast, feat.correlation, feat.energy, feat.homogeneity), ...
            'Interpreter','none');
    end
end
