% Elle Stark January 2024

% Script to calculate integral length scale and time scale for multisource
% plume dataset using function toolbox from:
% https://www.mathworks.com/matlabcentral/fileexchange/108944-a-note-on-the-calculation-of-the-length-scales-of-turbulence 
% Ref: Cheynet, E. (2016). Wind-induced vibrations of a suspension bridge: A case study in full-scale. PhD thesis. University of Stavanger, Norway.


clc; clear vars; close all;

% Set up H5 directory folders
h5_dir = 'D:/';
h5_file = [h5_dir 'singlesource_2d_extended/Re100_0_5mm_50Hz_singlesource_2d.h5'];

%% Read in numeric grids and velocity arrays from h5 file
startFrame = 1;
endFrame = 9001;

dt = h5read(h5_file,'/Model Metadata/timeResolution'); dt = 1/dt; % time step of model data, 50 Hz for COMSOL sims packaged in h5 structures
spatialResolution = h5read(h5_file,'/Model Metadata/spatialResolution');

xMesh = h5read(h5_file,'/Model Metadata/xGrid'); x_gridVector = min(xMesh,[],'all'):spatialResolution:max(xMesh,[],'all');
yMesh = h5read(h5_file,'/Model Metadata/yGrid'); y_gridVector = min(yMesh,[],'all'):spatialResolution:max(yMesh,[],'all');

% uStack = h5read(h5_file,'/Flow Data/u',[1 1 startFrame],[size(xMesh,1) size(xMesh,2) endFrame]); % About 10 GB RAM per variable (for all time)
dataStack = h5read(h5_file,'/Flow Data/u',[1 1001 startFrame],[size(xMesh,1) 500 endFrame]);
% vStack = h5read(h5_file,'/Flow Data/v',[1 1 startFrame],[size(xMesh,1) size(xMesh,2) endFrame]);

% odorStack = h5read(h5_file,'/Odor Data/c1a',[1 1 startFrame],[size(xMesh,1) size(xMesh,2) endFrame]);

timeArray = h5read(h5_file,'/Model Metadata/timeArray');

% Assign desired field to dataStack variable
% dataStack = uStack;

% Flatten spatial arrays into stack of vectors to input to functions
data_vecs = reshape(dataStack, [], size(dataStack, 3), 1);
[Ny, N] = size(data_vecs);
data_mean = mean(data_vecs, 2);

%% Compute integral length scale (L) in the streamwise direction

% L_test = Lx(data_vecs(100, :), data_mean(100), dt, 'method', 1);

% Preallocate arrays for Lx at each point
Lux = zeros(1, Ny);
Tux = zeros(1, Ny);

% Compute length scale with function in Lx.m using autocovariance function.
% Method 1 computes area from zero lag until first zero crossing; Method 2 uses an exponential fit.
for i=1:Ny
    [Lux(i), Tux(i)] = Lx(data_vecs(i, :), data_mean(i), dt, 'method', 1);
    if i/1000==0
        pct = i / Ny * 100;
        fprintf('%p percent complete. \n', pct);
    end
end

% reshape stack of vectors into 2D arrays
% Lux_array = reshape(Lux, size(dataStack, 1), size(dataStack, 2));
Tux_array = reshape(Tux, size(dataStack, 1), size(dataStack, 2));

%% Save and plot results

% Save arrays of integral time scales
save('Plots/Lux_array_u_extendedsim3.mat', 'Tux_array')
% save('Plots/Tux_array_v.mat', 'Tux_array')

% Plot integral length scales across domain
% figure;
% pcolor(xMesh, yMesh, Lux_array);
% shading flat;
% colormap('jet');
% colorbar;
% % caxis([0, 0.18]);
% title('Integral Length Scale (m) of v field', 'method: temporal autocorrelation * mean velocity');

% Plot integral time scales across domain
figure; 
pcolor(xMesh(:, 501:1000), yMesh(:, 501:1000), Tux_array);
shading interp;
axis equal;
axis tight;
set(gca, 'TickDir', 'out');
% colormap('Blues');
c = colorbar;
c.TickDirection = 'out';
lims = get(c, 'Limits');
c.Ticks = 0:0.25:2.5;
box off;
% c.Ticks = linspace(0, 0.7, 8);
% set(c, 'YTick', [0:0.1:0.7]);
caxis([0, 2.5]);
title('Integral Time Scale (s) of u field');
print(gcf, 'Plots/ITS_u3.png', '-dpng', '-r600');

