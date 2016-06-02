%Output results in VTK format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%Last modified: 25 may 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteADMResiduals2VTK(Directory, Problem, iter, Grid, Residual, Residualc, CoarseGrid, maxLevel)
%Write a VTK file
fileID = fopen(strcat(Directory,'/VTK/',Problem,'Residual',num2str(iter),'.vtk'), 'w');
fprintf(fileID, '# vtk DataFile Version 2.0\n');
fprintf(fileID, strcat(Problem, ' results: Matteo Simulator\n'));
fprintf(fileID, 'ASCII\n');
fprintf(fileID, '\n');
fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx +1, Grid.Ny+1, 2);
fprintf(fileID, '\n');
fprintf(fileID, ['X_COORDINATES ' num2str(Grid.Nx+1) ' float\n']);
fprintf(fileID, '%f ', 0:Grid.dx:Grid.Lx);
fprintf(fileID, '\n');
fprintf(fileID, ['Y_COORDINATES ' num2str(Grid.Ny+1) ' float\n']);
fprintf(fileID, '%f ', 0:Grid.dy:Grid.Ly);
fprintf(fileID, '\n');
fprintf(fileID, 'Z_COORDINATES 2 float\n');
fprintf(fileID, '%d ', [0 1]);
fprintf(fileID, '\n');
fprintf(fileID, '\n');
fprintf(fileID, 'CELL_DATA   %d\n', Grid.N);
fprintf(fileID, '\n');
%FineScale residual
PrintScalar2VTK(fileID, Residual(1:Grid.N), ' FineResidualNw');
fprintf(fileID, '\n');
PrintScalar2VTK(fileID, Residual(Grid.N+1:end), ' FineResidualW');
fprintf(fileID, '\n');
%ADM residual
PrintScalar2VTK(fileID, Residualc(1:Grid.N), ' CoarseResidualNw');
fprintf(fileID, '\n');
PrintScalar2VTK(fileID, Residualc(Grid.N+1:end), ' CoarseResidualW');
fprintf(fileID, '\n');
%ADD ADM coarse grids
    PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
    for i=1:maxLevel
        fclose(fileID);
        fileID = fopen(strcat(Directory,'/VTK/',Problem,'Residual',num2str(i),'Level',num2str(iter),'.vtk'), 'w');
        fprintf(fileID, '# vtk DataFile Version 2.0\n');
        fprintf(fileID, strcat(Problem, ' results: Matteo Simulator\n'));
        fprintf(fileID, 'ASCII\n');
        fprintf(fileID, '\n');
        fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
        fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', CoarseGrid(i).Nx +1, CoarseGrid(i).Ny+1, 2);
        fprintf(fileID, '\n');
        fprintf(fileID, ['X_COORDINATES ' num2str(CoarseGrid(i).Nx+1) ' float\n']);
        fprintf(fileID, '%f ', 0:Grid.dx*3^i:Grid.Lx);
        fprintf(fileID, '\n');
        fprintf(fileID, ['Y_COORDINATES ' num2str(CoarseGrid(i).Ny+1) ' float\n']);
        fprintf(fileID, '%f ', 0:Grid.dy*3^i:Grid.Ly);
        fprintf(fileID, '\n');
        fprintf(fileID, 'Z_COORDINATES 2 float\n');
        fprintf(fileID, '%d ', [0 1]);
        fprintf(fileID, '\n');
        fprintf(fileID, '\n');
        fprintf(fileID, 'CELL_DATA   %d\n', CoarseGrid(i).Nx*CoarseGrid(i).Ny);
        PrintScalar2VTK(fileID, CoarseGrid(i).Active, ' ActiveCoarse');
        fprintf(fileID, '\n');
    end
fclose(fileID);    
end


