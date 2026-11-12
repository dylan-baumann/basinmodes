function [Z, X, Y, header] = read_esri_asc(filename, flip_y)
% READ_ESRI_ASC Parses an ESRI ASCII Raster file (.asc / .txt)
%
% Syntax:
%   [Z, X, Y, header] = read_esri_asc('bathymetry.asc')
%
% Inputs:
%   filename - Path to the ESRI ASCII raster file
%   flip_y   - (Boolean) whether or not to flip the returned raster so as to not be upside down.
%
% Outputs:
%   Z        - Matrix of bathymetry/elevation values (size: nrows x ncols)
%   X, Y     - Meshgrid matrices or coordinate vectors for plotting
%   header   - Struct containing file metadata (ncols, nrows, xllcorner, etc.)

    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    % Clean up file handle on function exit
    cleanup = onCleanup(@() fclose(fid));

    header = struct();
    
    % Parse the 6-line standard header
    for i = 1:6
        line = fgetl(fid);
        tokens = strsplit(strtrim(line));
        key = lower(tokens{1});
        val = str2double(tokens{2});
        
        switch key
            case 'ncols'
                header.ncols = val;
            case 'nrows'
                header.nrows = val;
            case {'xllcorner', 'xllcenter'}
                header.xll = val;
                header.x_type = key;
            case {'yllcorner', 'yllcenter'}
                header.yll = val;
                header.y_type = key;
            case 'cellsize'
                header.cellsize = val;
            case 'nodata_value'
                header.nodata = val;
            otherwise
                warning('Unknown header parameter: %s', key);
        end
    end
    
    % Read data matrix
    % ESRI ASCII files store data row-by-row from Top to Bottom
    Z = fscanf(fid, '%f', [header.ncols, header.nrows])';
    
    % Replace NO-DATA values with NaN for easy numerical processing
    if isfield(header, 'nodata')
        Z(Z == header.nodata) = NaN;
    end

    % Adjust origin if center coordinates were given instead of corner
    if strcmp(header.x_type, 'xllcenter')
        x0 = header.xll;
    else
        x0 = header.xll + (header.cellsize / 2);
    end

    if strcmp(header.y_type, 'yllcenter')
        y0 = header.yll;
    else
        y0 = header.yll + (header.cellsize / 2);
    end

    % Generate X and Y spatial coordinate vectors
    x_vec = x0 + (0:(header.ncols - 1)) * header.cellsize;
    
    % Y runs from top to bottom in the raster matrix
    y_vec = y0 + ((header.nrows - 1):-1:0) * header.cellsize;
    
    % Create spatial grids for plotting
    [X, Y] = meshgrid(x_vec, y_vec);

    if flip_y
        %X = flipud(X);  %shouldn't need, omit for perf reasons.
        Y = flipud(Y);
        Z = flipud(Z);
    end
end