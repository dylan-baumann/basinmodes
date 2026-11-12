function [Nv, VX, VY, K, EToV, BCType, node, edge] = mesh_polygon(poly, char_length, basename, gmsh_exe)
% MESH_POLYGON  Triangulate a straight-sided polygon with Gmsh and load it.
%
%   [Nv,VX,VY,K,EToV,BCType,node,edge] = mesh_polygon(poly, char_length, basename)
%   [...] = mesh_polygon(poly, char_length, basename, gmsh_exe)
%
%   poly         Nx2 polygon vertices (metres), e.g. from rect_geometry.m
%   char_length  Target triangle edge length (scalar, or one value per vertex)
%   basename     Stem for the generated files: writes <basename>.geo and
%                <basename>.msh in the current folder.
%   gmsh_exe     (optional) Gmsh command. Default 'gmsh'. Use a full path if
%                MATLAB/Octave cannot find it, e.g. '/opt/homebrew/bin/gmsh'.
%
%   The outputs are exactly the variables that readmsh.m returns, so they can
%   be dropped straight into basinmodes.m / basinmodes_curved.m. No .mat file
%   is needed: the mesh is built in memory, which is what we want when a
%   parameter sweep (e.g. peninsula size) rebuilds the mesh many times.
%
%   PIPELINE:   poly --polygon2geo--> .geo --Gmsh--> .msh --readmsh--> arrays
%
%   SANITY CHECKS (each one guards against a silent failure mode):
%     1. Mesh area equals polygon area (Gmsh meshed the right region).
%     2. Every boundary face is tagged Wall, and nothing else is. readmsh.m
%        tags boundaries geometrically via findedge(), which has caveats
%        (see the repo notes). A basin whose boundary faces are mis-tagged
%        would run without error but give wrong modes, so we cross-check the
%        tags against the mesh connectivity itself: a face is on the boundary
%        exactly when its neighbour element is itself (EToE(k,f) == k).
%        If the two disagree we warn and repair the tags from connectivity.

    WALL = 3;   % boundary-type code; must match Globals2D.m (Wall = 3)

    if nargin < 4 || isempty(gmsh_exe)
        gmsh_exe = 'gmsh';
    end

    geo_file = [basename '.geo'];
    msh_file = [basename '.msh'];

    % ---- 1. geometry -> .geo -------------------------------------------
    polygon2geo(poly, char_length, geo_file);

    % ---- 2. .geo -> .msh via Gmsh --------------------------------------
    % Delete any stale mesh so a failed Gmsh run can never silently load an
    % old file.
    if exist(msh_file, 'file')
        delete(msh_file);
    end
    cmd = sprintf('%s -2 -format msh22 -v 1 -o "%s" "%s"', gmsh_exe, msh_file, geo_file);
    [status, output] = system(cmd);
    if status ~= 0 || ~exist(msh_file, 'file')
        error(['mesh_polygon: Gmsh failed (exit status %d).\n' ...
               'Command: %s\n\nGmsh output:\n%s\n' ...
               'Is Gmsh installed and on the path seen by MATLAB/Octave? ' ...
               'Pass the full path as the 4th argument if not.'], ...
               status, cmd, output);
    end

    % ---- 3. .msh -> DG-FEM arrays --------------------------------------
    [Nv, VX, VY, K, EToV, BCType, node, edge] = readmsh(msh_file);

    % ---- 4. sanity checks ----------------------------------------------
    % (a) area: sum of triangle areas vs polygon area
    ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
    bx = VX(EToV(:,2)); by = VY(EToV(:,2));
    cx = VX(EToV(:,3)); cy = VY(EToV(:,3));
    tri_area  = 0.5 * abs((bx-ax).*(cy-ay) - (cx-ax).*(by-ay));
    mesh_area = sum(tri_area);
    poly_area = polyarea(poly(:,1), poly(:,2));
    rel_area_err = abs(mesh_area - poly_area) / poly_area;
    if rel_area_err > 1e-8
        warning('mesh_polygon: mesh area differs from polygon area by %.2e (relative).', rel_area_err);
    end

    % (b) boundary tagging vs connectivity
    [EToE, ~] = tiConnect2D(EToV);
    isBoundaryFace = (EToE == (1:K)' * ones(1, 3));   % neighbour is itself
    isTaggedWall   = (BCType == WALL);
    if any(isBoundaryFace(:) ~= isTaggedWall(:))
        warning(['mesh_polygon: readmsh tagged %d wall faces but the mesh has %d boundary faces. ' ...
                 'Re-tagging all boundary faces as Wall from connectivity.'], ...
                 nnz(isTaggedWall), nnz(isBoundaryFace));
        BCType = zeros(size(BCType));
        BCType(isBoundaryFace) = WALL;
    end

    % ---- 5. one-line summary -------------------------------------------
    Ea = [1 2 3]; Eb = [2 3 1];
    elen = sqrt((VX(EToV(:,Ea)) - VX(EToV(:,Eb))).^2 + ...
                (VY(EToV(:,Ea)) - VY(EToV(:,Eb))).^2);
    fprintf(['Mesh "%s": %d triangles, %d vertices, %d boundary faces; ' ...
             'edge length %.0f-%.0f m; area error %.1e\n'], ...
             basename, K, Nv, nnz(isBoundaryFace), min(elen(:)), max(elen(:)), rel_area_err);
end
