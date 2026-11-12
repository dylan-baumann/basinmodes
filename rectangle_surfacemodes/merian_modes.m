function [w, mn] = merian_modes(Lx, Ly, H, g, nmax)
% MERIAN_MODES  Analytic free-oscillation frequencies of a rectangular basin.
%
%   [w, mn] = merian_modes(Lx, Ly, H, g, nmax)
%
%   Closed rectangular basin (Lx by Ly), uniform depth H, NO rotation.
%   The surface displacement of mode (m,n) is
%
%       eta(x,y) ~ cos(m*pi*x/Lx) * cos(n*pi*y/Ly),     m,n = 0,1,2,...
%
%   (no-flux walls => zero normal derivative => cosines), with angular
%   frequency
%
%       omega_mn = pi * sqrt(g*H) * sqrt( (m/Lx)^2 + (n/Ly)^2 ).
%
%   For n = 0 this reduces to Merian's formula for a 1-D basin,
%   T_m = 2*Lx / (m*sqrt(g*H)).
%
%   Outputs (sorted by ascending frequency, (m,n) = (0,0) excluded):
%     w    column vector of angular frequencies  [rad/s]
%     mn   [m n] mode numbers for each entry of w
%   nmax   largest m and n to generate; it must be big enough that every
%          mode below the highest frequency you intend to compare is included.

    c = sqrt(g*H);
    [m, n] = meshgrid(0:nmax, 0:nmax);
    m = m(:);  n = n(:);
    keep = (m + n) > 0;                 % drop the trivial (0,0) "mode"
    m = m(keep);  n = n(keep);

    w = pi * c * sqrt((m/Lx).^2 + (n/Ly).^2);
    [w, order] = sort(w);
    mn = [m(order), n(order)];
end
