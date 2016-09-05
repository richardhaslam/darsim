%  ADM map 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 5 September 2016
%Last modified: 5 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_map < handle
    properties
       Nf % all cells belonging to l < lc
       Nc % all cells belonging to lc
       Nx % number of cells of level lc - 1 generated from lc
       Verteces 
       OriginalIndex
    end
    methods
        function Update(obj, FineGrid, CoarseGrid, ADMGrid)
            for level = 1:ADMGrid.MaxLevel
                
                % Number of cells
                obj.Nf(level) = sum(ADMGrid.N(1:level));
                obj.Nc(level) = ADMGrid.N(level + 1);
                obj.Nx(level) = ADMGrid.N(level + 1) * prod(CoarseGrid(level).CoarseFactor);
                
                % verteces 
                obj.Verteces = 1;
                
                % original indexes
                obj.OriginalIndex = 1;
            end
        end
    end
end