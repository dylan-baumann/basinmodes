function poly = rect_geometry(Lx, Ly)
% RECT_GEOMETRY  Corner points of an axis-aligned Lx-by-Ly rectangular basin.
%
%   poly = rect_geometry(Lx, Ly)
%
%   Returns a 4x2 array of [x y] vertices (metres), ordered counter-clockwise
%   starting from the south-west corner at the origin:
%
%        (0,Ly) o-----------o (Lx,Ly)
%               |           |
%        (0,0)  o-----------o (Lx,0)
%
%   The polygon is implicitly closed (last vertex connects back to the first).
%
%   This function is deliberately tiny. Every study site in this workflow is
%   "a function that returns a polygon". The next case (rectangle with a
%   peninsula) will be a sibling function that returns more vertices; nothing
%   downstream (polygon2geo -> Gmsh -> readmsh -> basinmodes) needs to change.

    if Lx <= 0 || Ly <= 0
        error('rect_geometry: Lx and Ly must be positive.');
    end

    poly = [0  0;
            Lx 0;
            Lx Ly;
            0  Ly];
end
