%Output results in VTK format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Write2VTK(Directory, Problem, timestep, Grid, K, P, S, CoarseGrid)
%Write a VTK file
fileID = fopen(strcat(Directory,'/VTK/',Problem,num2str(timestep - 1),'.vtk'), 'w');
fprintf(fileID, '# vtk DataFile Version 2.0\n');
fprintf(fileID, strcat(Problem, ' results: Matteo Simulator\n'));
fprintf(fileID, 'ASCII\n');
fprintf(fileID, '\n');
fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx +1, Grid.Ny+1, 2);
fprintf(fileID, '\n');
fprintf(fileID, ['X_COORDINATES ' num2str(Grid.Nx+1) ' float\n']);
fprintf(fileID, '%f ', 0:Grid.Nx);
fprintf(fileID, '\n');
fprintf(fileID, ['Y_COORDINATES ' num2str(Grid.Ny+1) ' float\n']);
fprintf(fileID, '%f ', 0:Grid.Ny);
fprintf(fileID, '\n');
fprintf(fileID, 'Z_COORDINATES 2 float\n');
fprintf(fileID, '%d ', [0 1]);
fprintf(fileID, '\n');
fprintf(fileID, '\n');
fprintf(fileID, 'CELL_DATA   %d\n', Grid.N);
fprintf(fileID, '\n');
%Permeability
PrintScalar2VTK(fileID, reshape(K(1,:,:), Grid.N, 1), ' PERMX');
fprintf(fileID, '\n');
%Pressure
PrintScalar2VTK(fileID, reshape(P, Grid.N, 1), ' PRESSURE');
fprintf(fileID, '\n');
%Saturation
PrintScalar2VTK(fileID, reshape(S, Grid.N, 1), ' SATURATION');
fprintf(fileID, '\n');
if (nargin > 7)
    %ADD ADM coarse grids
    PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
    fclose(fileID);
    fileID = fopen(strcat(Directory,'/VTK/',Problem,'Level1',num2str(timestep - 1),'.vtk'), 'w');
    fprintf(fileID, '# vtk DataFile Version 2.0\n');
    fprintf(fileID, strcat(Problem, ' results: Matteo Simulator\n'));
    fprintf(fileID, 'ASCII\n');
    fprintf(fileID, '\n');
    fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
    fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', CoarseGrid(1).Nx +1, CoarseGrid(1).Ny+1, 2);
    fprintf(fileID, '\n');
    fprintf(fileID, ['X_COORDINATES ' num2str(CoarseGrid(1).Nx+1) ' float\n']);
    fprintf(fileID, '%f ', 0:3:Grid.Nx);
    fprintf(fileID, '\n');
    fprintf(fileID, ['Y_COORDINATES ' num2str(CoarseGrid(1).Ny+1) ' float\n']);
    fprintf(fileID, '%f ', 0:3:Grid.Ny);
    fprintf(fileID, '\n');
    fprintf(fileID, 'Z_COORDINATES 2 float\n');
    fprintf(fileID, '%d ', [0 1]);
    fprintf(fileID, '\n');
    fprintf(fileID, '\n');
    fprintf(fileID, 'CELL_DATA   %d\n', CoarseGrid(1).Nx*CoarseGrid(1).Ny);
    PrintScalar2VTK(fileID, CoarseGrid(1).Active, ' ActiveCoarse');
    fprintf(fileID, '\n');
end
fclose(fileID);
end