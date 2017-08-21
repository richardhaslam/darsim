%  Prolongation builder MS pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 August 2017
%Last modified: 7 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef prolongation_builder_FAMSPressure < prolongation_builder_MSPressure
    properties
    end
    methods
        function obj = prolongation_builder_FAMSPressure(n, cf)
            obj@prolongation_builder_MSPressure(n, cf);
        end
        function ADMProlp = ADMProlongation(obj, ADMGrid, FineGrid, CoarseGrid, ADMRest)
            % Pressure prolongation
            ADMProlp = 1;
            
            % Loop over the levels
            for level = ADMGrid.MaxLevel:-1:2
               Prolp = obj.LevelProlongation(ADMGrid, CoarseGrid(level-1), level);
               % Multiply by previous objects
               ADMProlp = Prolp * ADMProlp;
            end
            
            % Last prolongation is different coz I use fine-scale ordering
            Prolp = obj.LastProlongation(ADMGrid, FineGrid);
            
            % Multiply by previous objects
            ADMProlp = Prolp * ADMProlp;
        end
        function Prolp = LevelProlongation(obj, ADMGrid, FineGrid, level)
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % For a given level the prolongation operator looks like this             
             %       Nf     Nc
             %    ----       ----
             %    |      |      |
             % Nf |  I   |  0   | 
             %    |______|______|           
             %    |      |      |
             % Nx |      |      |
             %    |      |      |
             %    ----       ----
             
             
             % Update map for the next level
             obj.ADMmap.Update(ADMGrid, FineGrid, level);
             
             % 1. Build the object
             Prolp = sparse(obj.ADMmap.Nf + obj.ADMmap.Nx, obj.ADMmap.Nf + obj.ADMmap.Nc);
             
             % 2. Fill in top left
             Prolp(1:obj.ADMmap.Nf, 1:obj.ADMmap.Nf) = speye(obj.ADMmap.Nf);
             
             % 3. Fill in Bottom left
             Prolp(obj.ADMmap.Nf + 1 : end,  obj.ADMmap.Verteces) = obj.P{level}(obj.ADMmap.OriginalIndexNx, obj.ADMmap.OriginalIndexVerteces);
             
             % 4. Fill in Bottom right
             Prolp(obj.ADMmap.Nf + 1 :end, obj.ADMmap.Nf + 1 : end) = obj.P{level}(obj.ADMmap.OriginalIndexNx, obj.ADMmap.OriginalIndexNc);
        end
        function Prolp = LastProlongation(obj, ADMGrid, FineGrid)
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % For the first level the prolongation operator looks like this             
             %       Nf     Nl1
             %    ----       ----
             %    |      |      |
             %    |      |   0  | 
             % Nl0|      |______|           
             %    |      |      |
             %    |      |      |
             %    |      |      |
             %    ----       ----
             
             % Update map
             obj.ADMmap.Update(ADMGrid, FineGrid, 1);
           
             % 1. Build object
             Prolp = sparse(FineGrid.N, obj.ADMmap.Nf + obj.ADMmap.Nc);
             
             % 2. Fill in FS verteces of level 1
             Prolp(:,  obj.ADMmap.Verteces) = obj.P{1}(:, obj.ADMmap.OriginalIndexVerteces);
             
             % 3. Fill in coarse-scale nodes
             Prolp(:, obj.ADMmap.Nf + 1 : end) = obj.P{1}(:, obj.ADMmap.OriginalIndexNc);
             
             % 4. Fill in fine-scale nodes first
             rows = obj.ADMmap.OriginalIndexNf';
             columns = 1:obj.ADMmap.Nf;
             Prolp(rows, :) = 0; % if it s fine-scale already I get rid of useless fillings
             Prolp(sub2ind(size(Prolp), rows, columns)) = 1;
        end
        function UpdateProlongationOperator(obj, FineGrid, CoarseGrid, ProductionSystem)
            % for now I do not update pressure basis functions
        end
    end
end