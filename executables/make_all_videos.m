function make_all_videos(caseName)

% Suppress all warnings
warning('off', 'all');

% Move from executables/ to Project/
cd ..

% Full path to case results
casePath = fullfile("postProcessing", caseName);

% Read domain info from 0.00/Domain.csv
domainFile = fullfile(casePath, "0.00", "Domain.csv");
domain = readmatrix(domainFile);
x_loc = domain(1, :);
y_loc = domain(2, :);
[X, Y] = meshgrid(x_loc, y_loc);

% Sort time folders
folderList = dir(casePath);
timeFolders = {};
for i = 1:length(folderList)
    name = folderList(i).name;
    if folderList(i).isdir && ~ismember(name, {'.', '..'})
        if ~isempty(regexp(name, '^\d+(\.\d+)?$', 'once'))
            timeFolders{end+1} = name;
        end
    end
end

% Sort numerically
timeValues = str2double(timeFolders);
[~, sortIdx] = sort(timeValues);
sortedFolders = timeFolders(sortIdx);

% Set frame rate
frameRate = 5;

% Create video writers
vu = VideoWriter(fullfile(casePath, ['u_velocity_' caseName '.mp4']), 'MPEG-4'); vu.FrameRate = frameRate;
vv = VideoWriter(fullfile(casePath, ['v_velocity_' caseName '.mp4']), 'MPEG-4'); vv.FrameRate = frameRate;
vs = VideoWriter(fullfile(casePath, ['streamfunction_' caseName '.mp4']), 'MPEG-4'); vs.FrameRate = frameRate;
vz = VideoWriter(fullfile(casePath, ['vorticity_' caseName '.mp4']), 'MPEG-4'); vz.FrameRate = frameRate;

open(vu); open(vv); open(vs); open(vz);

for i = 1:length(sortedFolders)
    folderName = sortedFolders{i};
    basePath = fullfile(casePath, folderName);

    try
        u = readmatrix(fullfile(basePath, "U.csv"));
        v = readmatrix(fullfile(basePath, "V.csv"));
        stream = readmatrix(fullfile(basePath, "StreamFunction.csv"));
        zeta = readmatrix(fullfile(basePath, "Vorticity.csv"));

        % --- u velocity ---
        figU = figure('Visible', 'off');
        contourf(X, Y, u, 100, 'LineColor', 'none');
        colorbar; clim([-1 1]); axis equal tight;
        writeVideo(vu, getframe(figU));
        close(figU);

        % --- v velocity ---
        figV = figure('Visible', 'off');
        contourf(X, Y, v, 100, 'LineColor', 'none');
        colorbar; clim([-1 1]); axis equal tight;
        writeVideo(vv, getframe(figV));
        close(figV);

        % --- stream function ---
        figS = figure('Visible', 'off');
        contourf(X, Y, stream, 100, 'LineColor', 'none');
        colorbar; axis equal tight;
        writeVideo(vs, getframe(figS));
        close(figS);

        % --- vorticity ---
        figZ = figure('Visible', 'off');
        contourf(X, Y, zeta, 100, 'LineColor', 'none');
        colorbar; axis equal tight;
        writeVideo(vz, getframe(figZ));
        close(figZ);

    catch
        % Silent skip on missing/corrupt data
    end
end

% Close videos
close(vu); close(vv); close(vs); close(vz);

% Reset warnings
warning('on', 'all');
end
