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
            fileID = fopen(strcat(obj.FileName, num2str(obj.VTKindex),'.vtk'), 'w');
            fprintf(fileID, '# vtk DataFile Version 2.0\n');
            fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
            if obj.isBinary
                fprintf(fileID, 'BINARY\n');
            else
            	fprintf(fileID, 'ASCII\n');
            end
            fprintf(fileID, '\n');
            
            fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');
            
            fprintf(fileID, ['POINTS ' num2str(Grid.N*8) ' double\n']);
            Corners = horzcat(Grid.CornerPointGridData.Cell.SW_Bot_Corner,...
                              Grid.CornerPointGridData.Cell.SE_Bot_Corner,...
                              Grid.CornerPointGridData.Cell.NW_Bot_Corner,...
                              Grid.CornerPointGridData.Cell.NE_Bot_Corner,...
                              Grid.CornerPointGridData.Cell.SW_Top_Corner,...
                              Grid.CornerPointGridData.Cell.SE_Top_Corner,...
                              Grid.CornerPointGridData.Cell.NW_Top_Corner,...
                              Grid.CornerPointGridData.Cell.NE_Top_Corner)';
            Corners = reshape(Corners(:),3,Grid.N*8)';
            Corners(:,3) = -Corners(:,3);
            fprintf(fileID, '%f %f %f\n',Corners');
            fprintf(fileID, '\n');
            
            fprintf(fileID, ['CELLS ' num2str(Grid.N) ' ' num2str(Grid.N*(1+8)) '\n']);
            indMatrix = horzcat(8*ones(Grid.N,1) , reshape(0:Grid.N*8-1, 8,Grid.N)' );
            fprintf(fileID, '%d %d %d %d %d %d %d %d %d\n',indMatrix');
            fprintf(fileID, '\n');
            
            fprintf(fileID, ['CELL_TYPES ' num2str(Grid.N) '\n']);
            fprintf(fileID, '%d ', 11*ones(1,Grid.N));
            fprintf(fileID, '\n\n');
            
            % Print all existing variables
            fprintf(fileID, 'CELL_DATA   %d\n', Grid.N);
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
            
            %ADD ADM ACTIVETime
            obj.PrintScalar2VTK(fileID, Grid.ActiveTime, ' ACTIVETime');
            fprintf(fileID, '\n');
            
            %ADD ADM ACTIVEFine (coarse grids)
            obj.PrintScalar2VTK(fileID, Grid.Active, ' ACTIVEFine');
            fprintf(fileID, '\n');
            
            fclose(fileID);
        end
    end
    methods (Access = private)
        function PrintScalar2VTK(obj, fileID, scalar, name)
            %Print a scalar in VTK format
            fprintf(fileID, ' \n');
            fprintf(fileID, strcat('SCALARS  ', name,' float 1\n'));
            fprintf(fileID, 'LOOKUP_TABLE default\n');
            if obj.isBinary
                fwrite(fileID, scalar','float', 'b');
            else
            	fprintf(fileID,'%d ', scalar);
            end
        end
        function PrintVector2VTK(obj, fileID, vector, name)
            %Print a vector in VTK format
            fprintf(fileID, ' \n');
            fprintf(fileID, strcat('VECTORS  ', name,' float \n'));
            if obj.isBinary
                fwrite(fileID, vector','float', 'b');
            else
            	fprintf(fileID,'%d ', vector);
            end
        end
    end
end