function [xxr, yyr, rst_resampled] = resample_raster(xx, yy, rst, resample_resolution)
    xmin = min(xx(:)); xmax = max(xx(:)); xlen = xmax - xmin;
    ymin = min(yy(:)); ymax = max(yy(:)); ylen = ymax - ymin;

    dx = resample_resolution; dy = dx;
    Nx = ceil(xlen/dx); Ny = ceil(ylen/dy);
    xr = xmin + (0:Nx-1)*dx;
    yr = ymin + (0:Ny-1)*dy;
    [xxr, yyr] = meshgrid(xr, yr);

    rst_resampled = interp2(xx, yy, rst, xxr, yyr, 'linear', 0);
end