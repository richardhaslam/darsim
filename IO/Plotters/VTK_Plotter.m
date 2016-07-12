% VTK Plotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 11 July 2016
%Last modified: 11 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef VTK_Plotter < Plotter
    properties
        FileName
        VTKindex 
    end
    methods
        function obj = VTK_Plotter(Directory, Problem)
            obj.FileName = strcat(Directory, '/Output/VTK/', Problem);
            obj.VTKindex = 0;
        end
        function PlotWells(obj, Inj, Prod)
            %InJectors
            for i=1:length(Inj)
                name = strcat(obj.FileName,num2str(i),'Inj','.vtk');
                WriteAWell(Inj(i), name, Grid);
            end 
            %Producers
            for i=1:length(Prod)
                name = strcat(obj.FileName, num2str(i),'Prod','.vtk');
                WriteAWell(Prod(i), name, Grid);
            end
            
        end
        function WriteAWell(Well, name, Grid)
            fileID = fopen(name, 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'Matteo Simulator: Well s file\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'ASCII\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET POLYDATA\n');
            fprintf(fileID, ['POINTS ', num2str(length(Well.cells) + 1), ' float\n']);
            %loop over perforations
            z1 = Grid.h/2;
            z2 = 2*Grid.h;
            for c=1:length(Well.cells)
                dummy = mod(Well.cells(c), Grid.Nx);
                if dummy == 0
                    i = Grid.Nx;
                else
                    i = dummy;
                end
                j = (Well.cells(c) - i)/Grid.Nx + 1;
                x = (i - 1)*Grid.dx + Grid.dx/2;
                y = (j - 1)*Grid.dy + Grid.dy/2;
                fprintf(fileID, '%10.5f %10.5f %10.5f\n', x, y, z1);
            end
            fprintf(fileID, '%10.5f %10.5f %10.5f\n', x, y, z2);
            fprintf(fileID, ['LINES ', num2str(length(Well.cells)), ' ',num2str(3*(length(Well.cells))), '\n']);
            for c = 1:length(Well.cells)
                fprintf(fileID, '%d %d %d\n', [2, c-1, c]);
            end
            fclose(fileID);
        end
        function PlotSolution(obj, Status, Grid)
            %Write a VTK file
            fileID = fopen(strcat(obj.FileName, num2str(obj.VTKindex),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
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
            %ADD ADM coarse grids
            PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
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
            obj.VTKindex = obj.VTKindex + 1;
        end
        function PlotPermeability(obj, Grid, K)
            %Permeability
            fileID = fopen(strcat(obj.FileName, num2str(0),'.vtk'), 'a');
            PrintScalar2VTK(fileID, reshape(K(1,:,:), Grid.N, 1), ' PERMX');
            fprintf(fileID, '\n');
        end
        function PlotResidual(obj)
        end
        function PlotBasisFunctions(obj)
        end
        function PlotADMGrid(obj, Grid, CoarseGrid)
            for i=1:maxLevel
                fileID = fopen(strcat(obj.FileName, num2str(i),'Level', num2str(obj.VTKindex),'.vtk'), 'w');
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
                fclose(fileID);
            end
        end
    end
    methods (Access = private)
        function PrintScalar2VTK(fileID, scalar, name)
            %Print a scalar in VTK format
            fprintf(fileID, strcat('SCALARS ', name,' float 1\n'));
            fprintf(fileID, 'LOOKUP_TABLE default\n');
            fprintf(fileID,'%d ', scalar);
        end
    end
end

%% Notes: 
% I have to change the spacing for coarse grids. CoarseGrid should be aware
% of the deltaX and deltaY of each grid.