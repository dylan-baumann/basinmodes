
%read-in gmesh .msh mesh format. mesh should already be triangulated
%using gmesh. 

%spits out mesh data and BC data for use in a NUDG code.
%note: you may need to tweak which BC flag. i.e. wall, neumann, dirichlet etc.

%filename = 'opeongonw.msh';

function [Nv,VX, VY, K, EToV,BCType,node,edge] = readmsh(filename)

    fid = fopen(filename,'r');

    %skip over meshformat header and nodes header
    fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);

    numNodes = str2num(fgetl(fid));

    Nv= numNodes;

    for jj=1:numNodes
        dump = str2num(fgetl(fid));
        VX(jj) = dump(2);
        VY(jj) = dump(3);
    end

    
    fgetl(fid); % '$EndNodes
    fgetl(fid); % '$Elements''

    numElements = str2num(fgetl(fid)); 
    K = numElements;

    node = [];

    %read list of boundary nodes
    readingnodes = true;
    while readingnodes == true
        dump = str2num(fgetl(fid));

        if dump(2) ~= 15  %we hit an entry that's not a node.
            readingnodes = false;

            edge(1) = dump(end-1);
            edge(2) = dump(end);
        else  %we're reading boundary nodes
            mynode = [VX(dump(end)) VY(dump(end))];
            node = [node; mynode];
        end
    end  %done reading nodes

    readingedges = true;
    edges=0;
    while readingedges == true
        dump = str2num(fgetl(fid));

        if dump(2) ~= 1 %we hit an entry that's not an edge
            readingedges = false;

            EToV(1) = dump(end-2);
            EToV(2) = dump(end-1);
            EToV(3) = dump(end);
        else  %we're reading edges
            myedge(1) = dump(end-1);
            myedge(2) = dump(end);

            edge = [edge; myedge];

        end
    end

    readingelements = true;
    while readingelements == true
        dumpstr = fgetl(fid);

        if strcmp(dumpstr,'$EndElements') == 1 %end of element list
            readingelements = false;
        else %we're reading elements
            dump = str2num(dumpstr);

            myEToV(1) = dump(end-2);
            myEToV(2) = dump(end-1);
            myEToV(3) = dump(end);

            EToV = [EToV; myEToV];
        end

    end

    K = length(EToV);

    %we've got all the info we need, so get out of dodge.
    fclose(fid);


    % Reorder elements to ensure counter clockwise orientation
    ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
    bx = VX(EToV(:,2)); by = VY(EToV(:,2));
    cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

    D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
    i = find(D<0);
    EToV(i,:) = EToV(i,[1 3 2]);
    %done reordering

    % Build connectivity matrix
    [EToE, EToF] = tiConnect2D(EToV);

    Vert = [VX' VY'];

    %allocate BCType table.
    BCType = 0*EToF;
    
    
    
    %the 'node' data used by gmsh refers to its own spline data
    %not the vertices of the finite element mesh.

    node=[VX' VY'];
    

    %Insert the correct BC codes for boundaries
    Wall=3;
    BCType = CorrectBCTable_v3(EToV,VX,VY,node,edge,BCType,Wall, 1e-2);
end