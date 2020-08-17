% VTK Plotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef CornerPointGrid_VTK_Plotter < VTK_Plotter
    properties
    end
    methods
        function WriteAWell(obj, Well, name, Grid)
        end
        function PlotReservoirSolution(obj, Reservoir, Grid)
            %Write a VTK file for Reservoir
            fileID = fopen(strcat(obj.FileName, '_Reservoir', num2str(obj.VTKindex,'%04d'),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            if obj.isBinary
                fprintf(fileID, 'BINARY\n');
            else
            	fprintf(fileID, 'ASCII\n');
            end
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');
            
            fprintf(fileID, ['POINTS ' num2str(Grid.CornerPointGridData.N_Nodes) ' double\n']);
            fprintf(fileID, '%f %f %f\n',Grid.CornerPointGridData.Nodes');
            fprintf(fileID, '\n');
            
            fprintf(fileID, ['CELLS ' num2str(Grid.N) ' ' num2str(Grid.N*(1+8)) '\n']);
            indMatrix = horzcat(8*ones(Grid.N,1) , Grid.CornerPointGridData.Cells.Vertices-1 );
            fprintf(fileID, '%d %d %d %d %d %d %d %d %d\n',indMatrix');
            fprintf(fileID, '\n');
            
            fprintf(fileID, ['CELL_TYPES ' num2str(Grid.N) '\n']);
            fprintf(fileID, '%d ', 11*ones(1,Grid.N));
            fprintf(fileID, '\n\n');
            
            % Print all existing variables
            fprintf(fileID, 'CELL_DATA %d\n', Grid.N);
            N_var = double(Reservoir.State.Properties.Count);
            Names = Reservoir.State.Properties.keys;
            for i=1:N_var
                if strcmp(Reservoir.State.Properties(Names{i}).Type, 'scalar')
                    obj.PrintScalar2VTK(fileID, Reservoir.State.Properties(Names{i}).Value, [' ',Names{i}]);
                else
                    obj.PrintVector2VTK(fileID, Reservoir.State.Properties(Names{i}).Value, [' ',Names{i}]);
                end
                fprintf(fileID, '\n');
            end

            % Add the "fractured" flag for reservoir grid cells that are overlapped by a fracture (if any)
            FracturedFlag = zeros(Grid.N,1);
            if ~isempty(Grid.ListOfFracturedReservoirCells)
                FracturedFlag(Grid.ListOfFracturedReservoirCells) = 1;
            end
            obj.PrintScalar2VTK(fileID, FracturedFlag, ' IsFractured');
            fprintf(fileID, '\n');
            
            % Add the cell volume
            obj.PrintScalar2VTK(fileID, Grid.Volume, ' Volume');

            % Add ADM ACTIVETime
            obj.PrintScalar2VTK(fileID, Grid.ActiveTime, ' ACTIVETime');
            fprintf(fileID, '\n');
            
            % Add ADM ACTIVEFine (coarse grids)
            obj.PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
            fprintf(fileID, '\n');
            
            % Add manual flag for R5 model
            if Grid.N == 1276945
                Flagged = zeros(Grid.N,1);
                Cell_Indices = [959588;1029937;1033734;1193170;1266538];
                Flagged(Cell_Indices)=1;
                for n = 1 : length(Cell_Indices)
                    Flagged( Grid.CornerPointGridData.Cell.Index_Neighbors{n} ) = 2;
                end
                obj.PrintScalar2VTK(fileID, Flagged, ' isFlagged');
                fprintf(fileID, '\n');
            end
            
            fclose(fileID);
            
            obj.PlotInternalFaces(Reservoir, Grid);
        end
        function PlotInternalFaces(obj, Reservoir, Grid)
            %Write a VTK file for Reservoir
            fileID = fopen(strcat(obj.FileName, '_Interfaces', num2str(obj.VTKindex,'%04d'),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            if obj.isBinary
                fprintf(fileID, 'BINARY\n');
            else
            	fprintf(fileID, 'ASCII\n');
            end
            fprintf(fileID, '\n');
            fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');
            
            fprintf(fileID, ['POINTS ' num2str(Grid.CornerPointGridData.N_Nodes) ' double\n']);
            fprintf(fileID, '%f %f %f\n',Grid.CornerPointGridData.Nodes');
            fprintf(fileID, '\n');
            
            N_InternalFaces = Grid.CornerPointGridData.N_InternalFaces;
            maxNumVertices = max(Grid.CornerPointGridData.Internal_Faces.N_Vertices);
            FaceVertices = cell2mat( cellfun(@(x) [x ones(1,maxNumVertices-numel(x))*x(end)] , Grid.CornerPointGridData.Internal_Faces.Vertices , 'uni' , 0 ) );
            fprintf(fileID, ['CELLS ' num2str(N_InternalFaces) ' ' num2str(N_InternalFaces*(1+maxNumVertices)) '\n']);
            
            indMatrix = horzcat(maxNumVertices*ones(N_InternalFaces,1) , FaceVertices-1 );
            
            FormatSpec =[];
            for n = 1 : maxNumVertices
                FormatSpec = strcat(FormatSpec,"%d ");
            end
            FormatSpec = char(FormatSpec);
            FormatSpec(end)=[];
            FormatSpec = strcat(FormatSpec,'\n');
            fprintf(fileID, FormatSpec,indMatrix');
            fprintf(fileID, '\n');
            
            fprintf(fileID, ['CELL_TYPES ' num2str(N_InternalFaces) '\n']);
            fprintf(fileID, '%d ', 9*ones(1,N_InternalFaces));
            fprintf(fileID, '\n\n');
            
            % Print all existing variables
            fprintf(fileID, 'CELL_DATA %d\n', N_InternalFaces);
            N_var = 1;
            Names{1} = 'Trans';
            for i=1:N_var
                obj.PrintScalar2VTK(fileID, Grid.Trans, [' ',Names{i}]);
                fprintf(fileID, '\n');
            end

            fclose(fileID);
        end
    end
    methods (Access = private)
        function PrintScalar2VTK(obj, fileID, scalar, name)
            %Print a scalar in VTK format
            fprintf(fileID, ' \n');
            fprintf(fileID, strcat('SCALARS  ', name,' double 1\n'));
            fprintf(fileID, 'LOOKUP_TABLE default\n');
            if obj.isBinary
                fwrite(fileID, scalar','double', 'b');
            else
            	fprintf(fileID,'%1.5e ', scalar);
            end
        end
        function PrintVector2VTK(obj, fileID, vector, name)
            %Print a vector in VTK format
            fprintf(fileID, ' \n');
            fprintf(fileID, strcat('VECTORS  ', name,' double \n'));
            if obj.isBinary
                fwrite(fileID, vector','double', 'b');
            else
            	fprintf(fileID,'%1.5e ', vector);
            end
        end
    end
end