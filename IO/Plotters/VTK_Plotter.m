% VTK Plotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 11 July 2016
%Last modified: 17 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef VTK_Plotter < Plotter
    properties
        FileName
        VTKindex 
    end
    methods
        function obj = VTK_Plotter(Directory, Problem)
            if ~exist(strcat(Directory,'/Output/VTK/'), 'dir')
                mkdir(Directory,'/Output/VTK');
            else
                delete(strcat(Directory,'/Output/VTK/*.vtk'));
            end
            obj.FileName = strcat(Directory, '/Output/VTK/', Problem);
            obj.VTKindex = 0;
        end
        function PlotWells(obj, Inj, Prod, Grid)
            %InJectors
            for i=1:length(Inj)
                name = strcat(obj.FileName,num2str(i),'Inj','.vtk');
                obj.WriteAWell(Inj(i), name, Grid);
            end 
            %Producers
            for i=1:length(Prod)
                name = strcat(obj.FileName, num2str(i),'Prod','.vtk');
                obj.WriteAWell(Prod(i), name, Grid);
            end
            
        end
        function WriteAWell(obj, Well, name, Grid)
            fileID = fopen(name, 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'Matteo Simulator: Well s file\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'ASCII\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET POLYDATA\n');
            fprintf(fileID, ['POINTS ', num2str(length(Well.Cells) + 1), ' float\n']);
            %loop over perforations
            z1 = Grid.dz/2;
            z2 = 2*Grid.dz;
            for c=1:length(Well.Cells)
                dummy = mod(Well.Cells(c), Grid.Nx);
                if dummy == 0
                    i = Grid.Nx;
                else
                    i = dummy;
                end
                j = (Well.Cells(c) - i)/Grid.Nx + 1;
                x = (i - 1)*Grid.dx + Grid.dx/2;
                y = (j - 1)*Grid.dy + Grid.dy/2;
                fprintf(fileID, '%10.5f %10.5f %10.5f\n', x, y, z1);
            end
            fprintf(fileID, '%10.5f %10.5f %10.5f\n', x, y, z2);
            fprintf(fileID, ['LINES ', num2str(length(Well.Cells)), ' ',num2str(3*(length(Well.Cells))), '\n']);
            for c = 1:length(Well.Cells)
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
            fprintf(fileID, '%f ', 0:Grid.dx:Grid.dx * Grid.Nx);
            fprintf(fileID, '\n');
            fprintf(fileID, ['Y_COORDINATES ' num2str(Grid.Ny+1) ' float\n']);
            fprintf(fileID, '%f ', 0:Grid.dy:Grid.dy * Grid.Ny);
            fprintf(fileID, '\n');
            fprintf(fileID, 'Z_COORDINATES 2 float\n');
            fprintf(fileID, '%d ', [0 1]);
            fprintf(fileID, '\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'CELL_DATA   %d\n', Grid.N);
            fprintf(fileID, '\n');
            %ADD ADM coarse grids
            obj.PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
            fprintf(fileID, '\n');
            %Pressure
            obj.PrintScalar2VTK(fileID, Status.p, ' PRESSURE');
            fprintf(fileID, '\n');
            %Saturation
            obj.PrintScalar2VTK(fileID, Status.S, ' SATURATION');
            fprintf(fileID, '\n');
            %Capillary Pressure
            obj.PrintScalar2VTK(fileID, Status.pc, ' CapPRESSURE');
            fprintf(fileID, '\n');
            %x11
            obj.PrintScalar2VTK(fileID, Status.x1(:,1), ' x1w');
            fprintf(fileID, '\n');
            %x12
            obj.PrintScalar2VTK(fileID, Status.x1(:,2), ' x1nw');
            fprintf(fileID, '\n');
            %x21
            obj.PrintScalar2VTK(fileID, 1-Status.x1(:,1), ' x2w');
            fprintf(fileID, '\n');
            %x22
            obj.PrintScalar2VTK(fileID, 1-Status.x1(:,2), ' x2nw');
            fprintf(fileID, '\n');
            %z1
            obj.PrintScalar2VTK(fileID, Status.z(:,1), ' z1');
            fprintf(fileID, '\n');
            %Density
            obj.PrintScalar2VTK(fileID, Status.rho(:,1), ' rhoW');
            fprintf(fileID, '\n');
            obj.PrintScalar2VTK(fileID, Status.rho(:,2), ' rhoNw');
            fprintf(fileID, '\n');
            obj.PrintScalar2VTK(fileID, Status.rho(:,1).*Status.S + Status.rho(:,2).*(1 - Status.S), ' rhoT');
            obj.VTKindex = obj.VTKindex + 1;
        end
        function PlotPermeability(obj, Grid, K)
            %Permeability
            fileID = fopen(strcat(obj.FileName, num2str(0),'.vtk'), 'a');
            obj.PrintScalar2VTK(fileID, reshape(K(:,1), Grid.N, 1), ' PERMX');
            fprintf(fileID, '\n');
        end
        function PlotResidual(obj)
        end
        function PlotBasisFunctions(obj)
        end
        function PlotADMGrid(obj, Grid, CoarseGrid)
            for i=1:length(CoarseGrid)
                fileID = fopen(strcat(obj.FileName, num2str(i),'Level', num2str(obj.VTKindex),'.vtk'), 'w');
                fprintf(fileID, '# vtk DataFile Version 2.0\n');
                fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
                fprintf(fileID, 'ASCII\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
                fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', CoarseGrid(i).Nx +1, CoarseGrid(i).Ny+1, 2);
                fprintf(fileID, '\n');
                fprintf(fileID, ['X_COORDINATES ' num2str(CoarseGrid(i).Nx+1) ' float\n']);
                fprintf(fileID, '%f ', 0 : Grid.dx * CoarseGrid(i).CoarseFactor(1) : Grid.dx * Grid.Nx);
                fprintf(fileID, '\n');
                fprintf(fileID, ['Y_COORDINATES ' num2str(CoarseGrid(i).Ny+1) ' float\n']);
                fprintf(fileID, '%f ', 0 : Grid.dy * CoarseGrid(i).CoarseFactor(2) : Grid.dy * Grid.Ny);
                fprintf(fileID, '\n');
                fprintf(fileID, 'Z_COORDINATES 2 float\n');
                fprintf(fileID, '%d ', [0 1]);
                fprintf(fileID, '\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'CELL_DATA   %d\n', CoarseGrid(i).N);
                obj.PrintScalar2VTK(fileID, CoarseGrid(i).Active, ' ActiveCoarse');
                fprintf(fileID, '\n');
                fclose(fileID);
            end
        end
    end
    methods (Access = private)
        function PrintScalar2VTK(obj, fileID, scalar, name)
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