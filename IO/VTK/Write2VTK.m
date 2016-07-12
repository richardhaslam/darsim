%Output results in VTK format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last modified: 30 june 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Write2VTK(Directory, Problem, timestep, Grid, K, Status, ADMActive, CoarseGrid, maxLevel, basisfunction)
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
PrintScalar2VTK(fileID, Status.p, ' PRESSURE');
fprintf(fileID, '\n');
%Saturation
PrintScalar2VTK(fileID, Status.s, ' SATURATION');
fprintf(fileID, '\n');
%Capillary Pressure
PrintScalar2VTK(fileID, Status.pc, ' CapPRESSURE');
fprintf(fileID, '\n');
%x11
PrintScalar2VTK(fileID, Status.x1(:,1), ' x1w');
fprintf(fileID, '\n');
%x12
PrintScalar2VTK(fileID, Status.x1(:,2), ' x1nw');
fprintf(fileID, '\n');
%x21
PrintScalar2VTK(fileID, 1-Status.x1(:,1), ' x2w');
fprintf(fileID, '\n');
%x22
PrintScalar2VTK(fileID, 1-Status.x1(:,2), ' x2nw');
fprintf(fileID, '\n');
%z1
PrintScalar2VTK(fileID, Status.z(:,1), ' z1');
%Density
PrintScalar2VTK(fileID, Status.rho(:,1), ' rhoW');
PrintScalar2VTK(fileID, Status.rho(:,2), ' rhoNw');
PrintScalar2VTK(fileID, Status.rho(:,1).*Status.s + Status.rho(:,2).*(1 - Status.s), ' rhoT');

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