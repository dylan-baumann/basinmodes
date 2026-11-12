% rectangle_surfacemodes.m
%
% Free surface oscillations (seiche modes) of a closed rectangular basin of
% uniform depth, computed with the basinmodes DG-FEM solver and checked
% against the analytic (Merian) solution.
%
% HOW TO RUN (from MATLAB or GNU Octave):
%     >> cd rectangle_surfacemodes
%     >> rectangle_surfacemodes
% The script must be run from its own folder (like the other case drivers)
% because it locates the repository root with cd('../').
%
% REQUIREMENT: Gmsh must be installed (command-line 'gmsh'). See GMSH_EXE below.
%
% WHAT THIS SCRIPT DOES, IN ORDER
%   1. Parameters      all physical / mesh / numerical choices in one place
%   2. Geometry+mesh   rect_geometry -> polygon2geo -> Gmsh -> readmsh
%   3. Solve           basinmodes.m (the shared engine; not edited)
%   4. Pick modes      separate gravity ("seiche") modes from vorticity modes
%   5. Validate        compare periods with the analytic Merian solution
%   6. Plot            first six seiche modes
%
% WHY basinmodes.m AND NOT basinmodes_curved.m
%   The curved engine fits a periodic spline through the boundary vertices and
%   bends the boundary elements onto it. For a smooth lake that is what you
%   want; for a rectangle it would round off the four corners. Straight-sided
%   triangles represent a rectangle EXACTLY, so the straight engine is both
%   simpler and more accurate here.

clear; close all;
currDir = pwd;
cd ../
root = pwd;
addpath(root);
setPaths;
cd(currDir);

Globals2D;   % must come BEFORE the mesh is loaded: the DG code shares
             % VX, VY, EToV, BCType, N, K, ... as global variables.

%% ======================== 1. PARAMETERS ================================

% --- Geometry (metres) -------------------------------------------------
% Balaton-like scale, but deliberately NOT Balaton's 10:1 shape: a squatter
% rectangle puts genuinely two-dimensional modes (m>=1 AND n>=1) among the
% first twenty, which makes a better test of the 2-D solver.
Lx = 60e3;          % length (x-direction)
Ly = 16e3;          % width  (y-direction)
H0 = 3.2;           % uniform depth (m)

% --- Physics -----------------------------------------------------------
g = 9.81;
LATITUDE = 46.8;    % degrees north (used only if ROTATION = true)
ROTATION = true;   % false: Coriolis effectively off -> exact Merian test
                    % true : Coriolis at LATITUDE (no closed-form answer here)
if ROTATION
    % Same expression the other drivers use.
    f = 2*(2*pi/(3600*24))*sin(LATITUDE*pi/180);
else
    % Negligible rather than exactly zero, so that sigma/f (used in the
    % titles of scrollmodes.m / animate_mode.m) is large instead of 0/0.
    f = 1e-9;
end

% --- Mesh --------------------------------------------------------------
MESH_CHARACTERISTIC_LENGTH = 2000;   % target triangle edge length (m)
GMSH_EXE = 'gmsh';                
                                    

% --- Numerics ----------------------------------------------------------
N = 2;                 % DG polynomial order inside each triangle
numpot  = 60;          % # potential (gravity-mode) basis functions kept
numstrm = 60;          % # streamfunction (vorticity-mode) basis functions kept
NUM_H_REFINES = 0;     % uniform mesh refinements applied inside the engine
NCOMPARE = 20;         % how many modes to compare against the analytic solution
TOL_PASS = 1e-2;       % pass if max relative period error is below this

% --- Flags the engine expects -----------------------------------------
USEMEANDEPTH   = false;
ANALYTIC_DEPTH = true;
H_analytic     = @(x, y) H0 + 0*x;    % uniform depth (must accept arrays)
DUMP_TO_FILE   = false;               % true -> engine saves basinmodes_output.mat

%% ======================== 2. GEOMETRY AND MESH =========================
poly = rect_geometry(Lx, Ly);
[Nv, VX, VY, K, EToV, BCType, node, edge] = ...
    mesh_polygon(poly, MESH_CHARACTERISTIC_LENGTH, 'rectangle', GMSH_EXE);

fprintf('DG degrees of freedom per field: %d (N=%d, K=%d)\n', K*(N+1)*(N+2)/2, N, K);

%% ======================== 3. SOLVE =====================================
% Everything the engine needs is now in the workspace. It builds the DG
% operators, solves the two Laplacian eigenproblems (Neumann -> potential,
% Dirichlet -> streamfunction), assembles the coupled rotating problem and
% solves it. See basinmodes.m. (The long list of numbers it prints is the
% vector of periods in hours: it is missing a semicolon in the engine.)
basinmodes;

fprintf('\nBasin area  : %.6e m^2  (Lx*Ly = %.6e)\n', A, Lx*Ly);
fprintf('Mean depth  : %.6f m     (H0      = %.6f)\n', Hbar, H0);

%% ======================== 4. PICK THE SEICHE MODES =====================
% The engine returns 2*numpot + numstrm eigenvalues in 'myeigs', sorted by
% absolute frequency. Their meaning:
%   * the first numstrm entries   : slow rotation-dominated "vorticity" modes
%                                   (|omega| ~ f or below; ~0 if f ~ 0)
%   * all remaining entries       : surface gravity modes = the seiches.
%                                   Each physical mode appears TWICE, as +omega
%                                   and -omega (complex-conjugate fields).
% So to get one entry per physical seiche we keep the entries with positive
% frequency. 'seiche_idx' is the index of each such mode in myeigs / eta, in
% order of increasing frequency: this vector is the "mode picker".

abs_w = abs(myeigs(:));
if max(abs_w(1:numstrm)) >= min(abs_w(numstrm+1:end))
    warning(['Vorticity and gravity bands overlap: the "skip the first numstrm ' ...
             'modes" rule used by the repo scripts is unreliable for this run.']);
end

tail = (numstrm+1):numel(myeigs);
seiche_idx = tail(real(myeigs(tail)) > 0);
if numel(seiche_idx) ~= numpot
    warning('Expected %d positive-frequency gravity modes, found %d.', numpot, numel(seiche_idx));
end

omega = real(myeigs(seiche_idx));     % rad/s, ascending
T_h   = 2*pi ./ omega / 3600;         % periods in hours

%% ======================== 5. VALIDATE vs MERIAN ========================
% Two independent checks, because right eigenvalues do not by themselves
% guarantee right eigenfunctions (and the spatial structure is the point):
%   (a) period of each mode vs the analytic Merian period
%   (b) shape of each mode vs the analytic cos(m*pi*x/Lx)*cos(n*pi*y/Ly):
%       relative L2 misfit after fitting a free amplitude (and sign).
%       (Integrals use nodal products, so this is a diagnostic accurate to
%       roughly the discretisation error, not an exact norm.)
% Both assume the k-th numerical mode corresponds to the k-th analytic mode,
% which holds when no two modes have nearly the same frequency (true here).
[w_mer, mn_mer] = merian_modes(Lx, Ly, H0, g, 60);
T_mer_h = 2*pi ./ w_mer / 3600;

nc = min([NCOMPARE, numel(T_h), numel(T_mer_h)]);
rel_err = (T_h(1:nc) - T_mer_h(1:nc)) ./ T_mer_h(1:nc);

shape_err = nan(nc, 1);
if ~ROTATION
    for k = 1:nc
        e = mode_snapshot(eta{seiche_idx(k)});
        a = cos(mn_mer(k,1)*pi*x/Lx) .* cos(mn_mer(k,2)*pi*y/Ly);
        s = dgint(e.*a, V, J) / dgint(a.*a, V, J);          % best-fit amplitude
        shape_err(k) = sqrt(dgint((e - s*a).^2, V, J) / dgint(e.^2, V, J));
    end
end

fprintf('\n   #   (m,n)   T_Merian [h]   T_numeric [h]   rel. error   shape misfit\n');
for k = 1:nc
    fprintf('  %2d   (%d,%d)   %11.5f   %12.5f   %+10.2e   %10.2e\n', ...
            k, mn_mer(k,1), mn_mer(k,2), T_mer_h(k), T_h(k), rel_err(k), shape_err(k));
end

if ~ROTATION
    fprintf('\nMax |relative period error| over first %d modes: %.2e\n', nc, max(abs(rel_err)));
    fprintf('Max mode-shape misfit      over first %d modes: %.2e\n', nc, max(shape_err));
    if max(abs(rel_err)) < TOL_PASS && max(shape_err) < TOL_PASS
        fprintf('PASS (tolerance %.0e)\n', TOL_PASS);
    else
        fprintf('FAIL (tolerance %.0e): check mesh size, N, numpot, or geometry.\n', TOL_PASS);
    end
else
    fprintf(['\nROTATION = true: Merian is no longer exact (f = %.3e 1/s), so the ' ...
             '"error" column is the shift in period caused by rotation, and the ' ...
             'shape column is NaN (no analytic shapes).\n'], f);
end

%% ======================== 6. PLOT ======================================
% mode_snapshot rotates the arbitrary eigenvector phase so each panel shows
% the standing wave at maximum displacement (plain real(eta) can be blank).
% Coordinates are passed in km purely for readable axes.
figure(2); clf; colormap(darkjet);
nplot = min(6, numel(seiche_idx));
for k = 1:nplot
    subplot(3, 2, k);
    PlotField2D(N, x/1e3, y/1e3, mode_snapshot(eta{seiche_idx(k)}));
    view([0 90]); axis equal; axis tight; colorbar;
    if ~ROTATION && k <= size(mn_mer, 1)
        title(sprintf('(m,n)=(%d,%d)   T=%.3f h', mn_mer(k,1), mn_mer(k,2), T_h(k)));
    else
        title(sprintf('mode %d   T=%.3f h', k, T_h(k)));
    end
    if k >= nplot-1, xlabel('distance (km)'); end
    drawnow;
end

%% ======================== RESULTS STRUCT ===============================
% Everything a later sweep script needs, in one place.
results.T_h        = T_h;
results.omega      = omega;
results.seiche_idx = seiche_idx;
results.K          = K;
results.rel_err    = rel_err;
results.shape_err  = shape_err;
