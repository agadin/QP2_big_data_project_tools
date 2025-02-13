image_analysis_and_mass_estimation
function image_analysis_and_mass_estimation
% IMAGE_ANALYSIS_AND_MASS_ESTIMATION
% This MATLAB script performs two tasks:
% 1. It scans a user‐selected directory for subfolders whose names contain 
%    either “CT” or “MRI” (as specified by the user). In each such folder it 
%    processes each *.tif image to compute:
%       - the image center (of non‐black pixels) [mm],
%       - minimum and maximum radii (distance from center) [mm],
%       - standard deviation of nonblack pixel intensities,
%       - and area (in mm²).
%    The results are saved to a CSV file.
%
% 2. It also computes, for each image, an estimated mass (with propagated 
%    error) using a pixel‐by‐pixel density (linearly interpolated between 
%    density_dark and density_bright). Imaging parameters (pixel spacing and 
%    slice thickness) are used to compute a voxel volume (and its uncertainty).
%    The per‐image mass (and error) are saved to a separate CSV file, and the
%    total mass per folder as well as an overall average mass are printed.
%
% (Adapt the density values below to suit your data.)

    %% --- User Input & Directory Selection ---
    selected_dir = uigetdir(pwd, 'Select Directory containing CT and MRI folders');
    if isequal(selected_dir, 0)
        disp('No directory selected.');
        return;
    end

    scan_type = upper(strtrim(input('Enter scan type to process (CT or MRI): ', 's')));
    if ~ismember(scan_type, {'CT','MRI'})
        disp('Invalid input. Please enter "CT" or "MRI".');
        return;
    end

    %% --- PART 1: Image Analysis ---
    folderResults = []; % will hold a struct array with one entry per valid image
    dList = dir(selected_dir);
    for i = 1:length(dList)
        if dList(i).isdir && ~ismember(dList(i).name, {'.', '..'})
            if contains(upper(dList(i).name), scan_type)
                folder_path = fullfile(selected_dir, dList(i).name);
                fprintf('Processing folder: %s\n', dList(i).name);
                results = processFolder(folder_path);
                if ~isempty(results)
                    % Add the folder name to each record
                    for k = 1:length(results)
                        results(k).Folder = dList(i).name;
                    end
                    folderResults = [folderResults; results]; %#ok<AGROW>
                else
                    fprintf('  No valid images found in %s.\n', dList(i).name);
                end
            end
        end
    end

    if isempty(folderResults)
        disp('No valid images found in any target folders.');
    else
        output_csv = fullfile(selected_dir, 'image_analysis_output.csv');
        writeCSV(output_csv, folderResults);
    end

    %% --- PART 2: Mass Estimation ---
    % Set the density values (example values; adjust as needed):
    density_bright = 1;   % e.g. for the brightest pixels (g/cm³)
    density_dark   = 0.3; % e.g. for the darkest nonzero pixels (g/cm³)

    [massRecords, folderMassTotals, overallAvgMass, overallMassError] = ...
        massEstimation(selected_dir, scan_type, density_dark, density_bright);

    % Print total mass per folder
    fprintf('\nTotal Estimated Mass per Folder:\n');
    folderNames = fieldnames(folderMassTotals);
    for i = 1:length(folderNames)
        fprintf('  %s: %.3f g\n', folderNames{i}, folderMassTotals.(folderNames{i}));
    end
    fprintf('\nOverall Average Estimated Mass: %.3f g ± %.3f g\n', overallAvgMass, overallMassError);

    mass_csv = fullfile(selected_dir, 'mass_estimation_output.csv');
    massTable = struct2table(massRecords);
    writetable(massTable, mass_csv);
    fprintf('Mass estimation details saved to: %s\n', mass_csv);
end

%% ======================== LOCAL FUNCTIONS ================================

function result = processImage(image_path)
% PROCESSIMAGE loads a TIFF image, converts it to grayscale (if needed),
% and computes various statistics on non-black pixels.
%
% Returns a struct with fields:
%   - center:      [row_mm, col_mm] of the nonzero pixel centroid (mm)
%   - min_radius:  Minimum distance from the center (mm)
%   - max_radius:  Maximum distance from the center (mm)
%   - intensity_std: Standard deviation of nonzero intensities
%   - area:        Area (mm²) based on nonzero pixels

    PIXEL_SCALE = 10.0 / 17.53;  % mm per pixel

    % Try to load the image
    try
        data = imread(image_path);
        if ndims(data) == 3  % if RGB, convert to grayscale
            data = rgb2gray(data);
        end
        data = double(data);  % work in double precision
    catch ME
        fprintf('Warning: Failed to read %s. Error: %s\n', image_path, ME.message);
        result = [];
        return;
    end

    total_pixels = numel(data);
    black_pixels = sum(data(:) == 0);
    black_ratio = black_pixels / total_pixels;
    if black_ratio > 0.99
        % Skip if the image is almost entirely black.
        result = [];
        return;
    end

    % Find indices of nonzero pixels and compute the centroid.
    [rows, cols] = find(data ~= 0);
    center_row = mean(rows);
    center_col = mean(cols);
    center_pixels = [center_row, center_col];
    center_mm = center_pixels * PIXEL_SCALE;

    % Compute distances (in pixels) from the center and convert to mm.
    distances_pixels = sqrt((rows - center_row).^2 + (cols - center_col).^2);
    min_radius = min(distances_pixels) * PIXEL_SCALE;
    max_radius = max(distances_pixels) * PIXEL_SCALE;

    % Compute the standard deviation and area (in mm²) for non-black pixels.
    non_black_intensities = data(data ~= 0);
    intensity_std = std(non_black_intensities);
    area_pixels = numel(non_black_intensities);
    area_mm2 = area_pixels * (PIXEL_SCALE^2);

    % Store computed values in the output structure.
    result.center = center_mm;
    result.min_radius = min_radius;
    result.max_radius = max_radius;
    result.intensity_std = intensity_std;
    result.area = area_mm2;
end

function results = processFolder(folder_path)
% PROCESSFOLDER processes all *.tif images in the specified folder.
%
% It returns a struct array (one element per image) with image analysis
% results. Each struct will later have a 'Folder' field added.
    files = dir(fullfile(folder_path, '*.tif'));
    results = [];
    for idx = 1:length(files)
        image_path = fullfile(folder_path, files(idx).name);
        res = processImage(image_path);
        if ~isempty(res)
            % Save filename and image index (using zero‐based index as in Python)
            res.filename = files(idx).name;
            res.index = idx - 1;
            results = [results; res]; %#ok<AGROW>
        end
    end
end

function writeCSV(output_path, folderResults)
% WRITECSV writes the image analysis results to a CSV file.
%
% The CSV will have the following columns:
% Folder, Image Index, Filename, Center_Y (mm), Center_X (mm),
% Min_Radius (mm), Max_Radius (mm), Intensity_STD, Area (mm²)
    n = length(folderResults);
    Folder = cell(n,1);
    ImageIndex = zeros(n,1);
    Filename = cell(n,1);
    Center_Y = zeros(n,1);
    Center_X = zeros(n,1);
    Min_Radius = zeros(n,1);
    Max_Radius = zeros(n,1);
    Intensity_STD = zeros(n,1);
    Area = zeros(n,1);

    for i = 1:n
        Folder{i} = folderResults(i).Folder;
        ImageIndex(i) = folderResults(i).index;
        Filename{i} = folderResults(i).filename;
        % In our convention, center(1) is the row (Y) and center(2) is the column (X)
        Center_Y(i) = folderResults(i).center(1);
        Center_X(i) = folderResults(i).center(2);
        Min_Radius(i) = folderResults(i).min_radius;
        Max_Radius(i) = folderResults(i).max_radius;
        Intensity_STD(i) = folderResults(i).intensity_std;
        Area(i) = folderResults(i).area;
    end

    T = table(Folder, ImageIndex, Filename, Center_Y, Center_X, Min_Radius, Max_Radius, Intensity_STD, Area);
    writetable(T, output_path);
    fprintf('CSV file saved to: %s\n', output_path);
end

function [mass, mass_error] = computeImageMassAndError(image_path, density_dark, density_bright, pixel_spacing, slice_thickness, delta_p, delta_t, delta_density)
% COMPUTEIMAGEMASSANDERROR loads an image and computes the mass contribution of
% that slice (in grams) along with its propagated error.
%
% Parameters:
%   image_path      : Path to the TIFF image.
%   density_dark    : Density for the darkest nonzero pixels (g/cm³).
%   density_bright  : Density for the brightest pixels (g/cm³).
%   pixel_spacing   : Pixel spacing (mm per pixel).
%   slice_thickness : Slice thickness (mm).
%   delta_p         : Uncertainty in pixel spacing (mm).
%   delta_t         : Uncertainty in slice thickness (mm).
%   delta_density   : Uncertainty in density (g/cm³).
%
% Returns:
%   mass (in grams) and mass_error (in grams).

    try
        data = imread(image_path);
        if ndims(data) == 3
            data = rgb2gray(data);
        end
        data = double(data);
    catch ME
        fprintf('Error reading image %s: %s\n', image_path, ME.message);
        mass = 0.0;
        mass_error = 0.0;
        return;
    end

    mask = data > 0;
    if sum(mask(:)) == 0
        mass = 0.0;
        mass_error = 0.0;
        return;
    end

    % Map nonzero intensities (1-255) to a factor in [0,1]:
    factor = (data(mask) - 1) / 254.0;
    % Linear interpolation of density:
    pixel_densities = density_dark + factor * (density_bright - density_dark);

    % Compute voxel volume (convert from mm³ to cm³; 1 cm³ = 1000 mm³)
    voxel_volume_cm3 = (pixel_spacing^2 * slice_thickness) / 1000.0;
    mass = sum(pixel_densities(:)) * voxel_volume_cm3;

    % ----- Error Propagation -----
    % Relative error in voxel volume:
    rel_err_volume = sqrt((2 * delta_p / pixel_spacing)^2 + (delta_t / slice_thickness)^2);
    % Approximate relative error in effective density:
    avg_density = (density_dark + density_bright) / 2.0;
    rel_err_density = delta_density / avg_density;
    % Combine the relative errors:
    rel_err = sqrt(rel_err_volume^2 + rel_err_density^2);
    mass_error = mass * rel_err;
end

function [massRecords, folderMassTotals, overallAvgMass, overallMassError] = massEstimation(selected_dir, scan_type, density_dark, density_bright)
% MASSESTIMATION processes the target folders (those whose names contain
% scan_type) to compute per-image mass and error. It returns:
%   - massRecords: a struct array with one record per image.
%   - folderMassTotals: a struct with a field per folder holding the total mass.
%   - overallAvgMass: the average mass over all target folders.
%   - overallMassError: the overall mass uncertainty (combined in quadrature).

    % Parameters common to all images:
    pixel_spacing = 10.0 / 17.53;  % mm per pixel
    delta_p = 0.01;                % uncertainty in pixel spacing (mm)
    delta_density = 0.02;          % uncertainty in density (g/cm³)

    massRecords = [];
    total_mass_squared_error = 0.0;
    total_mass = 0.0;
    target_folder_count = 0;
    folderMassTotals = struct();

    dList = dir(selected_dir);
    for i = 1:length(dList)
        if dList(i).isdir && ~ismember(dList(i).name, {'.', '..'})
            if contains(upper(dList(i).name), scan_type)
                target_folder_count = target_folder_count + 1;
                folder_name = dList(i).name;
                folder_path = fullfile(selected_dir, folder_name);
                
                % Determine slice thickness and delta_t based on folder name:
                if contains(upper(folder_name), 'CT')
                    slice_thickness = 4.0;
                    delta_t = 0.2;
                elseif contains(upper(folder_name), 'MRI')
                    slice_thickness = 2.0;
                    delta_t = 0.1;
                else
                    slice_thickness = 1.0;
                    delta_t = 0.1;
                end

                files = dir(fullfile(folder_path, '*.tif'));
                folder_mass = 0.0;
                for j = 1:length(files)
                    image_path = fullfile(folder_path, files(j).name);
                    [mass, mass_error] = computeImageMassAndError(...
                        image_path, density_dark, density_bright, ...
                        pixel_spacing, slice_thickness, delta_p, delta_t, delta_density);
                    
                    rec.Folder = folder_name;
                    rec.ImageIndex = j - 1;
                    rec.Filename = files(j).name;
                    rec.SliceThickness = slice_thickness;
                    rec.Mass = mass;
                    rec.MassError = mass_error;
                    massRecords = [massRecords; rec]; %#ok<AGROW>
                    
                    folder_mass = folder_mass + mass;
                    total_mass_squared_error = total_mass_squared_error + mass_error^2;
                end
                folderMassTotals.(folder_name) = folder_mass;
                total_mass = total_mass + folder_mass;
            end
        end
    end

    if target_folder_count > 0
        overallAvgMass = total_mass / target_folder_count;
    else
        overallAvgMass = 0;
    end
    overallMassError = sqrt(total_mass_squared_error);
end
