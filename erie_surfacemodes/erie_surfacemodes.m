clear;
close all;
basinmodes_banner;
currDir = pwd;
cd ../
root = pwd;
addpath(root);
setPaths;
cd(currDir);

Globals2D;

if ~exist('ESRI_ASCII_INPUT_FILE', 'var')
    ESRI_ASCII_INPUT_FILE = 'lake_erie.asc';
end

if ~exist('MESH_CHARACTERISTIC_LENGTH', 'var')
    MESH_CHARACTERISTIC_LENGTH = 4000;
end

if ~exist('RESAMPLE_RESOLUTION', 'var')
    RESAMPLE_RESOLUTION = 1000;
end

disp('Pre-requisite: Checking for libgmsh installation...')

[status, output] = system('export LD_LIBRARY_PATH=/lib/x86_64-linux-gnu/; python -c "import ctypes; libgmsh = ctypes.CDLL(''libgmsh.so'')"');
if status ~= 0
    fprintf('A fatal exception has occurred:\n\n%s\n', output);
    fprintf('Try installing libgmsh-dev, or an equivalent package for your system.\n');
    return;
end
disp('libgmsh found.')

% Pre-processing:
% (1) load bathymetry from ESRI ASCII format.
flip_y = true;
disp('Parsing ASCII raster bathymetry file...')
[bath, xx, yy, h] = read_esri_asc(ESRI_ASCII_INPUT_FILE, flip_y);
disp('Done.')

bath_resampled = bath;
if RESAMPLE_RESOLUTION ~= h.cellsize
    fprintf('Resampling bathymetry raster to ''%d m''...\n', RESAMPLE_RESOLUTION);
    xmin = min(xx(:)); xmax = max(xx(:)); xlen = xmax - xmin;
    ymin = min(yy(:)); ymax = max(yy(:)); ylen = ymax - ymin;

    dx = RESAMPLE_RESOLUTION; dy = dx;
    Nx = ceil(xlen/dx); Ny = ceil(ylen/dy);
    xr = xmin + (0:Nx-1)*dx;
    yr = ymin + (0:Ny-1)*dy;
    [xxr, yyr] = meshgrid(xr, yr);

    bath_resampled = interp2(xx, yy, bath, xxr, yyr, 'linear', 0);
    disp('Done.');
end

% (2) Next, process the raster to vectorial boundary data using contour functionality.
disp('Extracting boundary contour shape data (raster2vector)...');
[node, edge, bdryspline, depthdata, edge_array] = bath2contour_spline(bath_resampled, RESAMPLE_RESOLUTION);
disp('Done.')

% (3) Output boundary data to "*.geo" file (Gmsh-native geometry format).
geo_file = strrep(ESRI_ASCII_INPUT_FILE, '.asc', '.geo');
disp('Serializing boundary geometry to GEO-file...');
bdry2geo(node, edge_array, MESH_CHARACTERISTIC_LENGTH, geo_file);
disp('Done.');

% (4) Have Gmsh process the geo-file into a 2D triangular mesh for us.
% ... Details: Call libgmsh-dev python wrapper script to get a .msh file (from the .geo file) in the right format.
disp('Meshing: Invoking Gmsh (Frontal-Delaunay; triangles)...')
system(sprintf([ ...
        'export LD_LIBRARY_PATH=/lib/x86_64-linux-gnu/; PYTHONPATH=%s /usr/bin/python3 -c "import geo2msh as g2m; g2m.generate_msh_from_geo(''%s'')"' ...
    ], root, geo_file ...
));
disp('Done.');

%------------------------------------------------------------------- 
% Now, msh-file exists.
% Main step -- read the mesh, and do the DG-FEM stuff.

msh_file = strrep(geo_file, '.geo', '.msh');
[Nv, VX, VY, K, EToV, BCType, node, edge] = readmsh(msh_file);

% Polynomial order used for approximation
N = 1;

K = length(EToV);
NUM_H_REFINES = 0;

USEMEANDEPTH = false; %use mean depth instead of full bathymetry?
ANALYTIC_DEPTH = false;
 
f = 2*(2*pi/(3600*24))*sin(42.15*pi/180); %43.7 deg. lattitude.
g = 9.81;

numpot = 100;  %number of potential basis functions
numstrm = 100; %number of streamfunction basis functions

DUMP_TO_FILE = true;

CURVED_IMPLEMENTATION = false;
%-------------------------------------------------------------------

if CURVED_IMPLEMENTATION == true
    basinmodes_curved;
else
    basinmodes;
end
