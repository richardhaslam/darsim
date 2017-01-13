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
            obj.VTKindex = 1;
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
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx +1, Grid.Ny+1, Grid.Nz+1);
            fprintf(fileID, '\n');
            fprintf(fileID, ['X_COORDINATES ' num2str(Grid.Nx+1) ' float\n']);
            fprintf(fileID, '%f ', 0:Grid.dx:Grid.dx * Grid.Nx);
            fprintf(fileID, '\n');
            fprintf(fileID, ['Y_COORDINATES ' num2str(Grid.Ny+1) ' float\n']);
            fprintf(fileID, '%f ', 0:Grid.dy:Grid.dy * Grid.Ny);
            fprintf(fileID, '\n');
            fprintf(fileID, ['Z_COORDINATES ' num2str(Grid.Nz+1) ' float\n']);
            fprintf(fileID, '%d ', 0:Grid.dz:Grid.dz * Grid.Nz);
            fprintf(fileID, '\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'CELL_DATA   %d\n', Grid.N);
            fprintf(fileID, '\n');
            %ADD ADM coarse grids
            obj.PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
            fprintf(fileID, '\n');
            %Pressure
            obj.PrintScalar2VTK(fileID, Status.p - Status.pc, ' Pph1');
            fprintf(fileID, '\n');
            obj.PrintScalar2VTK(fileID, Status.p, ' Pph2');
            fprintf(fileID, '\n');
            %Saturation
            Status.S(Status.S<1e-7) = 0;
            obj.PrintScalar2VTK(fileID, Status.S, ' Sph1');
            fprintf(fileID, '\n');
            %Capillary Pressure
            obj.PrintScalar2VTK(fileID, Status.pc, ' CapPRESSURE');
            fprintf(fileID, '\n');
            [~, nc] = size(Status.z);
            for c=1:nc
                %xcv
                obj.PrintScalar2VTK(fileID, Status.x(:,(c-1)*2+1), [' x',num2str(c),'ph1']);
                fprintf(fileID, '\n');
                %xcl
                obj.PrintScalar2VTK(fileID, Status.x(:,(c-1)*2+2), [' x',num2str(c),'ph2']);
                fprintf(fileID, '\n');
                %zc
                obj.PrintScalar2VTK(fileID, Status.z(:,c), [' z',num2str(c)]);
                fprintf(fileID, '\n');
            end
            %Density
            obj.PrintScalar2VTK(fileID, Status.rho(:,1), ' rhoPh1');
            fprintf(fileID, '\n');
            obj.PrintScalar2VTK(fileID, Status.rho(:,2), ' rhoPh2');
            fprintf(fileID, '\n');
            obj.PrintScalar2VTK(fileID, Status.rho(:,1).*Status.S + Status.rho(:,2).*(1 - Status.S), ' rhoT');
            obj.VTKindex = obj.VTKindex + 1;
        end
        function PlotPermeability(obj, Grid, K)
            %Permeability
            fileID = fopen(strcat(obj.FileName, num2str(obj.VTKindex - 1),'.vtk'), 'a');
            obj.PrintScalar2VTK(fileID, reshape(K(:,1), Grid.N, 1), ' PERMX');
            fprintf(fileID, '\n');
        end
        function PlotBasisFunctions(obj, Grid, CoarseGrid, Prolp)
            %% 1. Level 1
            fileID = fopen(strcat(obj.FileName,'_BF_Level1.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            fprintf(fileID, 'ASCII\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx +1, Grid.Ny+1, Grid.Nz+1);
            fprintf(fileID, '\n');
            fprintf(fileID, ['X_COORDINATES ' num2str(Grid.Nx+1) ' float\n']);
            fprintf(fileID, '%f ', 0:Grid.dx:Grid.dx * Grid.Nx);
            fprintf(fileID, '\n');
            fprintf(fileID, ['Y_COORDINATES ' num2str(Grid.Ny+1) ' float\n']);
            fprintf(fileID, '%f ', 0:Grid.dy:Grid.dy * Grid.Ny);
            fprintf(fileID, '\n');
            fprintf(fileID, ['Z_COORDINATES ' num2str(Grid.Nz+1) ' float\n']);
            fprintf(fileID, '%d ', 0:Grid.dz:Grid.dz * Grid.Nz);
            fprintf(fileID, '\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'CELL_DATA   %d\n', Grid.N);
            fprintf(fileID, '\n');
            for j = 1:CoarseGrid(1).N
                obj.PrintScalar2VTK(fileID, full(Prolp{1}(:,j)),strcat(' BF',num2str(j)));
                fprintf(fileID, '\n');
            end
            fclose(fileID);
            %% 2. Levels > 1
            for i=2:length(CoarseGrid)
                fileID = fopen(strcat(obj.FileName,'_BF_Level',num2str(i),'.vtk'), 'w');
                fprintf(fileID, '# vtk DataFile Version 2.0\n');
                fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
                fprintf(fileID, 'ASCII\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
                fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', CoarseGrid(i-1).Nx +1, CoarseGrid(i-1).Ny+1, CoarseGrid(i-1).Nz+1);
                fprintf(fileID, '\n');
                fprintf(fileID, ['X_COORDINATES ' num2str(CoarseGrid(i-1).Nx+1) ' float\n']);
                fprintf(fileID, '%f ', 0 : Grid.dx * CoarseGrid(i-1).CoarseFactor(1) : Grid.dx * Grid.Nx);
                fprintf(fileID, '\n');
                fprintf(fileID, ['Y_COORDINATES ' num2str(CoarseGrid(i-1).Ny+1) ' float\n']);
                fprintf(fileID, '%f ', 0 : Grid.dy * CoarseGrid(i-1).CoarseFactor(2) : Grid.dy * Grid.Ny);
                fprintf(fileID, '\n');
                fprintf(fileID, ['Z_COORDINATES ' num2str(CoarseGrid(i-1).Nz+1) ' float\n']);
                fprintf(fileID, '%f ', 0 : Grid.dz * CoarseGrid(i-1).CoarseFactor(3) : Grid.dz * Grid.Nz);
                fprintf(fileID, '\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'CELL_DATA   %d\n', CoarseGrid(i-1).N);
                for j = 1:CoarseGrid(i).N
                    obj.PrintScalar2VTK(fileID, full(Prolp{i}(:,j)), strcat(' BF',num2str(j)));
                    fprintf(fileID, '\n');
                end
                fclose(fileID);
            end
        end
        function PlotDynamicBasisFunctions(obj, Grid, Prolp)
            fileID = fopen(strcat(obj.FileName,'Dynamic_BF_', num2str(obj.VTKindex - 1),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            fprintf(fileID, 'ASCII\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx +1, Grid.Ny+1, Grid.Nz+1);
            fprintf(fileID, '\n');
            fprintf(fileID, ['X_COORDINATES ' num2str(Grid.Nx+1) ' float\n']);
            fprintf(fileID, '%f ', 0:Grid.dx:Grid.dx * Grid.Nx);
            fprintf(fileID, '\n');
            fprintf(fileID, ['Y_COORDINATES ' num2str(Grid.Ny+1) ' float\n']);
            fprintf(fileID, '%f ', 0:Grid.dy:Grid.dy * Grid.Ny);
            fprintf(fileID, '\n');
            fprintf(fileID, ['Z_COORDINATES ' num2str(Grid.Nz+1) ' float\n']);
            fprintf(fileID, '%d ', 0:Grid.dz:Grid.dz * Grid.Nz);
            fprintf(fileID, '\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'CELL_DATA   %d\n', Grid.N);
            fprintf(fileID, '\n');
            if obj.VTKindex > 2
                [~, col] = size(Prolp);
                for j = 1:col
                    obj.PrintScalar2VTK(fileID, full(Prolp(:,j)),strcat(' DynamicBF',num2str(j)));
                    fprintf(fileID, '\n');
                end
                fclose(fileID);
            end
        end
        function PlotADMGrid(obj, Grid, CoarseGrid)
            obj.VTKindex = obj.VTKindex - 1;
            for i=1:length(CoarseGrid)
                fileID = fopen(strcat(obj.FileName, num2str(i),'Level', num2str(obj.VTKindex),'.vtk'), 'w');
                fprintf(fileID, '# vtk DataFile Version 2.0\n');
                fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
                fprintf(fileID, 'ASCII\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
                fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', CoarseGrid(i).Nx +1, CoarseGrid(i).Ny+1, CoarseGrid(i).Nz+1);
                fprintf(fileID, '\n');
                fprintf(fileID, ['X_COORDINATES ' num2str(CoarseGrid(i).Nx+1) ' float\n']);
                fprintf(fileID, '%f ', 0 : Grid.dx * CoarseGrid(i).CoarseFactor(1) : Grid.dx * Grid.Nx);
                fprintf(fileID, '\n');
                fprintf(fileID, ['Y_COORDINATES ' num2str(CoarseGrid(i).Ny+1) ' float\n']);
                fprintf(fileID, '%f ', 0 : Grid.dy * CoarseGrid(i).CoarseFactor(2) : Grid.dy * Grid.Ny);
                fprintf(fileID, '\n');
                fprintf(fileID, ['Z_COORDINATES ' num2str(CoarseGrid(i).Nz+1) ' float\n']);
                fprintf(fileID, '%f ', 0 : Grid.dz * CoarseGrid(i).CoarseFactor(3) : Grid.dz * Grid.Nz);
                fprintf(fileID, '\n');
                fprintf(fileID, '\n');
                fprintf(fileID, 'CELL_DATA   %d\n', CoarseGrid(i).N);
                obj.PrintScalar2VTK(fileID, CoarseGrid(i).Active, ' ActiveCoarse');
                fprintf(fileID, '\n');
                fclose(fileID);
            end
            obj.VTKindex = obj.VTKindex + 1;
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