%Output results in VTK format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Write2VTK(Directory, Problem, timestep, Grid, K, P, S, Pc, ADMActive, CoarseGrid, maxLevel, basisfunction, Status)

P = reshape(Status.p,Grid.Nx,Grid.Ny);
S = reshape(Status.s,Grid.Nx,Grid.Ny);

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
%Permeability
PrintScalar2VTK(fileID, reshape(K(1,:,:), Grid.N, 1), ' PERMX');
fprintf(fileID, '\n');
%Pressure
PrintScalar2VTK(fileID, reshape(P, Grid.N, 1), ' PRESSURE');
fprintf(fileID, '\n');
%Saturation
PrintScalar2VTK(fileID, reshape(S, Grid.N, 1), ' SATURATION');
fprintf(fileID, '\n');
PrintScalar2VTK(fileID, reshape(Pc, Grid.N, 1), ' CapPRESSURE');
fprintf(fileID, '\n');
if (ADMActive == 1)
if (basisfunction == 1)
    Nc = CoarseGrid(1).Nx * CoarseGrid(1).Ny;
    for c = 1:Nc
        PrintScalar2VTK(fileID, full(CoarseGrid(1).MsP(:,c)), strcat(' BasisFunction', num2str(c)));
        fprintf(fileID, '\n');
    end
end
    %ADD ADM coarse grids
    PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
    for i=1:maxLevel
        fclose(fileID);
        fileID = fopen(strcat(Directory,'/VTK/',Problem,num2str(i),'Level',num2str(timestep - 1),'.vtk'), 'w');
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
        %Output Basis functions
        if (basisfunction == 1 && i < maxLevel)
            Nc = CoarseGrid(i+1).Nx * CoarseGrid(i+1).Ny;
            for c = 1:Nc
                PrintScalar2VTK(fileID, full(CoarseGrid(i+1).MsP(:,c)), strcat(' BasisFunction', num2str(c)));
                fprintf(fileID, '\n');
            end
        end
    end
end
fclose(fileID);
end