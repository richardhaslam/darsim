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
       OriginalIndexCoarse
       OriginalIndexFine
    end
    methods
        function Update(obj, FineGrid, CoarseGrid, ADMGrid, level)
                % Number of cells
                obj.Nf = sum(ADMGrid.N(1:level));
                obj.Nc = ADMGrid.N(level + 1);
                obj.Nx = ADMGrid.N(level + 1) * prod(CoarseGrid(level).CoarseFactor);
                
                % Verteces 
                obj.Verteces = find(ADMGrid.Verteces(:,level));
                
                % original indexes
                obj.OriginalIndexVerteces = ;
                obj.OriginalIndexNc = ;
                obj.OriginalIndexNf = ;
        end
    end
end