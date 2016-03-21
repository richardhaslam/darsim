%ADM - Build R and P for static coarse grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R, Pp, Ps, NewDLGRGrid] = BuildRandP(FineGrid, CoarseGrid, DLGRGrid, level, Nc, Nf, Nx)
% Construct R and P from level l+1 to level l

if (level ~= 1)
    %% Levels > 1
    R = [diag(ones(Nx,1),0), zeros(Nx,Nf-Nx); zeros(Nc-Nx,Nf)];
    NewDLGRGrid = UpdateDLGRGrid(DLGRGrid, FineGrid, Nc, Nf, Nx);
    start = tic;
    if CoarseGrid.constant ~= 1
        Pp = R';
        start1 = tic;
        %BF in columns of centers of the level
        centers = find(DLGRGrid.Centers(:,level));
        % Copy MsP in DLGR prolongation operator
        Pp(Nx+1:Nf,centers) = CoarseGrid.MsP(NewDLGRGrid.CellIndex(Nx+1:Nf), DLGRGrid.Father(centers, level));
        Pp(Nx+1:Nf, Nx+1:Nc) = CoarseGrid.MsP(NewDLGRGrid.CellIndex(Nx+1:Nf), DLGRGrid.CellIndex(Nx+1:Nc));
        part1 = toc(start1);
        start2 = tic;
        for i = Nx+1:Nc
            c2 = DLGRGrid.CellIndex(i);
            for j = Nx+1:Nf
                c1 = NewDLGRGrid.Father(j,level);
                if (c1 == c2)
                    %Restriction
                    R(i,j) = 1;
                end
            end
        end
        part2 = toc(start2);
    else
        %CONSTANT
        for i = Nx+1:Nc
            for j = Nx+1:Nf
                if (NewDLGRGrid.Father(j,level) == DLGRGrid.CellIndex(i))
                    %Restriction
                    R(i,j) = 1;
                end
            end
        end
        Pp = R';
    end
else
    %% Level 1    
    %Fine level (keeps original ordering)
    R = zeros(Nc, Nf);
    if CoarseGrid.constant ~= 1
    Pp = zeros(Nf, Nc);
    %Fine Nodes that are centres
    for r=1:Nx
        c = DLGRGrid.Father(r,1);
        if ( DLGRGrid.I(r)== CoarseGrid.I(c) && DLGRGrid.J(r) == CoarseGrid.J(c))
            Pp(:,r) = CoarseGrid.MsP(:,c);
        end
    end
    %Coarse Nodes
    for r=Nx+1:Nc
        Pp(:,r) = CoarseGrid.MsP(:,DLGRGrid.CellIndex(r));
    end
    for c=1:Nc
        if (DLGRGrid.level(c) == 0)
            %Active fine cells
            R(c, DLGRGrid.CellIndex(c)) = 1;
            %If fine a set the row to  0 0 0 0
            Pp(DLGRGrid.CellIndex(c),:) = zeros(Nc,1);
            Pp(DLGRGrid.CellIndex(c),c) = 1;
        else
            % Active coarse cells
            R(c,:) = CoarseGrid.MsR(DLGRGrid.CellIndex(c),:);
        end
    end
    else
        %CONSTANT
        for c=1:Nc
        if (DLGRGrid.level(c) == 0)
            %Active fine cells
            R(c, DLGRGrid.CellIndex(c)) = 1;
        else
            % Active coarse cells
            R(c,:) = CoarseGrid.MsR(DLGRGrid.CellIndex(c),:);
        end
        end
        Pp = R';
    end
end
%% Make them sparse and create Ps
    R = sparse(R);
    Pp = sparse(Pp);
    Ps = R';
end