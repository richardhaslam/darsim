% VTK Plotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 11 July 2016
%Last modified: 17 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Plot_VTK < handle
    properties
        FileName
        VTKindex 
    end
    methods
        function obj = Plot_VTK(Directory, Problem)
            if ~exist(strcat(Directory,'/VTK/'), 'dir')
                mkdir(Directory,'/VTK');
            else
                delete(strcat(Directory,'/VTK/*.vtk'));
            end
            obj.FileName = strcat(Directory, '/VTK/', Problem);
            obj.VTKindex = 1;
        end
        function PlotSolution(obj, Simulation)
            obj.PlotReservoirSolution(Simulation.Reservoir);
            for f = 1 : length(Simulation.Fractures)
                obj.PlotFracturesSolution(Simulation.Fractures(f), f);
            end
            obj.VTKindex = obj.VTKindex + 1;
        end
        function PlotReservoirSolution(obj, Reservoir)
            %Write a VTK file for Reservoir
            fileID = fopen(strcat(obj.FileName, '_Matrix_', num2str(obj.VTKindex),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            fprintf(fileID, 'ASCII\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Reservoir.NX+1, Reservoir.NY+1, Reservoir.NZ+1);
            fprintf(fileID, '\n');
            fprintf(fileID, ['X_COORDINATES ' num2str(Reservoir.NX+1) ' float\n']);
            fprintf(fileID, '%f ', 0:Reservoir.DX:Reservoir.DX * Reservoir.NX);
            fprintf(fileID, '\n');
            fprintf(fileID, ['Y_COORDINATES ' num2str(Reservoir.NY+1) ' float\n']);
            fprintf(fileID, '%f ', 0:Reservoir.DY:Reservoir.DY * Reservoir.NY);
            fprintf(fileID, '\n');
            fprintf(fileID, ['Z_COORDINATES ' num2str(Reservoir.NZ+1) ' float\n']);
            fprintf(fileID, '%d ', 0:Reservoir.DZ:Reservoir.DZ * Reservoir.NZ);
            fprintf(fileID, '\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'CELL_DATA   %d\n', Reservoir.NX*Reservoir.NY*Reservoir.NZ);
            fprintf(fileID, '\n');
            fclose(fileID);
        end
        function PlotFracturesSolution(obj, Fracture, f)
            %Write a VTK file for each
            fileID = fopen(strcat(obj.FileName, '_Fracture', num2str(f,'%02d'), '_', num2str(obj.VTKindex),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            fprintf(fileID, 'ASCII\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET STRUCTURED_GRID\n');
            fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Fracture.N_Length_AB+1, Fracture.N_Width_AD+1, 1);
            fprintf(fileID, '\n');
            fprintf(fileID, 'POINTS    %d   double\n', size(Fracture.GridCoords,1)*size(Fracture.GridCoords,2) );
            GridCoords_X = Fracture.GridCoords(:,:,1);
            GridCoords_Y = Fracture.GridCoords(:,:,2);
            GridCoords_Z = Fracture.GridCoords(:,:,3);
            fprintf(fileID, '%f %f %f\n' , [GridCoords_X(:),GridCoords_Y(:),GridCoords_Z(:)]' ); 
            fprintf(fileID, '\n');
            fprintf(fileID, '\n');
            fprintf(fileID, 'CELL_DATA %d\n', Fracture.N_Length_AB*Fracture.N_Width_AD);
            fprintf(fileID, '\n');
            fclose(fileID);
        end
    end
    methods (Access = private)
        function PrintScalar2VTK(obj, fileID, scalar, name)
            %Print a scalar in VTK format
            fprintf(fileID, strcat('SCALARS  ', name,' float 1\n'));
            fprintf(fileID, 'LOOKUP_TABLE default\n');
            fprintf(fileID,'%d ', scalar);
        end
    end
end