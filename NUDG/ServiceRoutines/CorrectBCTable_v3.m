function BCType = CorrectBCTable_v3(EToV, VX, VY, node, edge, BCType, BCcode, tol, hack)
    %Modified by D. Steinmoeller to not require the global vars.
    %and to work without a distance function

    %this script can properly address elements with nodes on more than one wall,
    %or two faces along a 'corner' boundary.


    if ~exist('hack','var')
        hack = false;
    end
    if ~exist('tol', 'var')
        tol = 1e-2;
    end
    BCType = zeros(size(BCType));

    pxc = 0.5*(VX(EToV)+VX(EToV(:,[2 3 1])));
    pyc = 0.5*(VY(EToV)+VY(EToV(:,[2 3 1])));

    midps = [pxc(:) pyc(:)];

    %try this hack to fix colinearity bug with 'findedge'...
    %put in dummy point:
    if hack == true
        midps = [midps;1e9 1e9];
        edgenum = findedge(midps,node,edge,tol); %normal code.
        edgenum = edgenum(1:end-1);
    % drop dummy point's edge/flag. end hack.
    else
        edgenum = findedge(midps,node,edge,tol);
    end
    %if the face's midpoint has an edgenumber, then it lies on the boundary!
    idx = edgenum ~=0;

    BCType(edgenum>0) = BCcode;
return
