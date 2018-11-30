%  ADM map 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_map < handle
    properties
       Nf % all cells belonging to l < lc
       Nc % all cells belonging to lc
       Nx % number of cells of level lc - 1 generated from lc
       CF
       Verteces 
       OriginalIndexVerteces
       OriginalIndexNc
       OriginalIndexNx
       OriginalIndexNf
    end
    methods
        function obj = adm_map(cf)
            obj.CF = cf;
        end
        function Update(obj, ADMGrid, FineGrid, level)
                % Number of cells
                obj.Nf = sum(ADMGrid.N(:,1:level), 1);
                obj.Nc = ADMGrid.N(:, level+1);
                obj.Nx = ADMGrid.N(:, level+1) .* obj.CF;
                % obj.Nf = sum(ADMGrid.N(1:level));
                % obj.Nc = ADMGrid.N(level+1);
                % obj.Nx = obj.Nc * obj.CF;
                
                % Original Indeces of Nf
                obj.OriginalIndexNf = ADMGrid.CellIndex(1:sum(obj.Nf));
                
                % Verteces 
                obj.Verteces = find(ADMGrid.Verteces(:,level));
                obj.OriginalIndexVerteces = ADMGrid.Fathers(obj.Verteces, level);
                
                % Coarse nodes that will be prolonged
                obj.OriginalIndexNc = ADMGrid.CellIndex(sum(obj.Nf)+1 : end);
                
                % Update ADMGrid
                if level > 1
                    ADMGrid.Update(obj.Nx, obj.Nf, obj.Nc, FineGrid);
                
                    % New Fine-Grid indeces
                    obj.OriginalIndexNx = ADMGrid.CellIndex(sum(obj.Nf)+1 : end);
                end
        end
    end
end