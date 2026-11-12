function ret = bdry2geo(node, edge_array, characteristic_length, outfile_name)
% Converts node and edge (node connectivity) boundary data
% into Gmsh .geo format and saves to 'outfile_name'.
% This will create Gmsh points, lines (splines), line loops and plane-surfaces.

    fid = fopen(outfile_name,'w+');
    clen = characteristic_length;

    numpoints = length(node);
    numedgeblocks = length(edge_array);

    zcoord = 0; %2D plane, zcoord is always zero

    % Write points.
    for jj = 1:numpoints
        fprintf(fid,'Point(%d) = {%0.3f, %0.3f, %0.3f, %0.3f};\n', jj, node(jj,1), node(jj,2), zcoord, clen);
    end

    spline_size=8; %5
    splinenum=1;
    lineloopinds = [];
    for ii = 1:numedgeblocks
        edgeblock = edge_array{ii};
        vertlist = [edgeblock(1:end,1); edgeblock(1,1)];

        %make last index in loop a multiple of spline-size
        %so we can adaptively choose the number of points in the last spline.
        lastj = floor(length(vertlist)/spline_size)*(spline_size)-spline_size;

        splineinds=[];
        for j=0:spline_size:lastj
            if j == lastj
                endind = length(vertlist);
            else
                endind = j+spline_size+1;
            end
            vertliststring = sprintf('%d,', vertlist(j+1:endind));
            vertliststring = vertliststring(1:end-1); %remove trailing comma

            fprintf(fid,'Spline(%d) = {%s};\n',splinenum,vertliststring);
            splineinds=[splineinds splinenum];

            splinenum=splinenum+1;
        end
        %done with that loop, so write a corresponding lineloop to file
        lineloopstring = sprintf('%d,', splineinds);
        lineloopstring = lineloopstring(1:end-1);
        fprintf(fid,'Line Loop(%d) = {%s};\n',splinenum,lineloopstring);
        lineloopinds = [lineloopinds splinenum];
        splinenum=splinenum+1;
    end

    % Now that all splines and loops are defined, can define a surface...
    planesurfstring = sprintf('%d,', lineloopinds);
    planesurfstring = planesurfstring(1:end-1);
    fprintf(fid,'Plane Surface(%d) = {%s};\n',splinenum,planesurfstring);
    splinenum=splinenum+1;
    %done.

    fclose(fid);
    ret = [];
end
