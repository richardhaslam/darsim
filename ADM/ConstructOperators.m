function [R, Pp, Ps] = ConstructOperators(FineGrid, CoarseGrid, DLGRGrid)
%Construct DLGR R and P
R = struct('matrix', {});
Pp = struct('matrix', {});
Ps = struct('matrix', {});
maxLevel = DLGRGrid.level(end);
%Other Levels
for i=maxLevel:-1:2
    [Nf, Nx] = NewNumberOfCells(DLGRGrid, i);
    [R(i).matrix, Pp(i).matrix, Ps(i).matrix, DLGRGrid] = BuildRandP(CoarseGrid(i-1), CoarseGrid(i), DLGRGrid, i, sum(DLGRGrid.N), Nf, Nx);    
end
%First level
[R(1).matrix, Pp(1).matrix, Ps(1).matrix] = BuildRandP(FineGrid, CoarseGrid(1), DLGRGrid, 1, sum(DLGRGrid.N), FineGrid.Nx*FineGrid.Ny, DLGRGrid.N(1));
end

function [Nf, Nx] = NewNumberOfCells(DLGRGrid, x)
%For now I use 9 that is CF between two levels
Nx = sum(DLGRGrid.N) - DLGRGrid.N(x+1);
Nf = Nx + DLGRGrid.N(x+1) * 9;
end