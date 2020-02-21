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
        CornerPointGridData
        Trans
        HeatTrans
    end
    methods
        function obj = corner_point_grid(ReservoirProperties)
            obj.CornerPointGridData = ReservoirProperties.CornerPointGridData;
            obj.N = ReservoirProperties.Grid.N_ActiveCells;
            obj.Active = ones(obj.N, 1);
            obj.ActiveTime = ones(obj.N, 1);
            obj.Trans = zeros(ReservoirProperties.CornerPointGridData.N_InternalFaces,1);
            obj.HeatTrans = zeros(ReservoirProperties.CornerPointGridData.N_InternalFaces,1);
        end
        function Initialize(obj, Reservoir)
            obj.ComputeRockTransmissibilities();
        end
        function ComputeRockTransmissibilities(obj)
            Perm = obj.CornerPointGridData.Permeability;
            CellNeighbor1Index = obj.CornerPointGridData.CellNeighbor1Index;
            CellNeighbor2Index = obj.CornerPointGridData.CellNeighbor2Index;
            CellNeighbor1Vec = obj.CornerPointGridData.CellNeighbor1Vec;
            CellNeighbor2Vec = obj.CornerPointGridData.CellNeighbor2Vec;
            Nvec = obj.CornerPointGridData.Nvec;
            Trans_Half_1 = sum( Perm(CellNeighbor1Index,1) .* CellNeighbor1Vec .* Nvec , 2 ) ./ sum( CellNeighbor1Vec .* CellNeighbor1Vec , 2 );
            Trans_Half_2 = sum( Perm(CellNeighbor2Index,1) .* CellNeighbor2Vec .* Nvec , 2 ) ./ sum( CellNeighbor2Vec .* CellNeighbor2Vec , 2 );
            obj.Trans = Trans_Half_1 .* Trans_Half_2 ./ (Trans_Half_1 + Trans_Half_2);
        end
        
    end
end