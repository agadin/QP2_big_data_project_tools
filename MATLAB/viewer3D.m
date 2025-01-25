folder_path = uigetdir(pwd, 'Select the folder containing the .tif file');

    if folder_path == 0
        disp('Folder selection canceled.');
    else
        disp(['Selected folder: ', folder_path]);
    end
imds = imageDatastore(folder_path, ...
    'FileExtensions', {'.tif'}, ...
    'IncludeSubfolders', false); 

if isempty(imds.Files)
    error('No .tif files found in the specified folder.');
end

info = imfinfo(imds.Files{1}); 
volumeData = zeros(info.Height, info.Width, numel(imds.Files), 'uint16'); 

k = 1;
while hasdata(imds)
    volumeData(:, :, k) = read(imds);
    k = k + 1;
end

if contains(folder_path, 'CT', 'IgnoreCase', true)
    zSpacing = 100; 
else
    zSpacing = 50; 
end
xyScaleFactor = 15.73 / 10; 

scaledVolume = imresize3(volumeData, ...
    [size(volumeData, 1) * xyScaleFactor, size(volumeData, 2) * xyScaleFactor, size(volumeData, 3) * (zSpacing / 10)]);

viewer = viewer3d();

volshow(scaledVolume, ...
    'Parent', viewer, ...               
    'Colormap', gray, ...               
    'RenderingStyle', 'VolumeRendering'); 

viewer.BackgroundColor = [0, 0, 0]; 
