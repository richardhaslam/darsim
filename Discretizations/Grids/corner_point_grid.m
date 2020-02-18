%  Cartesian grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author:
%TU Delft
%Created:
%Last modified:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef corner_point_grid < grid_darsim
    properties
        CellNeighbor1Vec
        CellNeighbor1Index
        CellNeighbor2Vec
        CellNeighbor2Index
        Nvec
        X_Y_Z
     
    end
    methods
        function obj = corner_point_grid(CellNeighbor1Vec, CellNeighbor2Vec, Nvec, X_Y_Z, CellNeighbor1Index, CellNeighbor2Index)     
            obj.hT1 = sum(CellNeighbor1Vec .* Nvec .* X_Y_Z(CellNeighbor1Index,:), 2)...
                ./ sum(CellNeighbor1Vec .* CellNeighbor1Vec, 2);
            
            obj.hT2 = abs(sum(CellNeighbor2Vec .* Nvec .* X_Y_Z(CellNeighbor2Index,:), 2)...
                ./ sum(CellNeighbor2Vec .* CellNeighbor2Vec, 2));
            
            obj.hTF = 1 ./ (1 ./ obj.hT1 + 1 ./ obj.hT2);
        end
    end
end