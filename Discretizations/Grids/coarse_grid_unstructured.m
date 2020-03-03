%  Coarse grid class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini & Mousa HosseiniMehr
%TU Delft
%Created: 26 July 2016
%Last modified: 27 March 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef coarse_grid_unstructured < coarse_grid
    properties
    end
    methods
        function obj = coarse_grid_unstructured() 
        end
        function BuildCoarseGrid(obj, FineGrid)            
            % This function will read the CoarseNodeIndices and creates the
            % coarse grids.
        end
        function AddGridCoordinates(obj, FineGrid)
        end
        function AssignNeighbours(obj)
        end
        function AssignNeighbours2D(obj)
        end
        function AssignNeighbours1D(obj)
        end
        function AddWells(obj, Inj, Prod)
        end
    end
end