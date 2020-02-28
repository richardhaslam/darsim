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
        Nx
        Ny
        Nz
    end
    methods
        function obj = corner_point_grid(ReservoirProperties)
            obj.CornerPointGridData = ReservoirProperties.CornerPointGridData;
            obj.Nx = ReservoirProperties.Grid.N(1);
            obj.Ny = ReservoirProperties.Grid.N(2);
            obj.Nz = ReservoirProperties.Grid.N(3);
            obj.N = ReservoirProperties.Grid.N_ActiveCells; % For CornerPointGrid, this is the number of cells actively involved, not dot_prod(Nx,Ny,Nz)
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
            CellNeighbor1Index = obj.CornerPointGridData.Internal_Face.CellNeighbor1Index;
            CellNeighbor2Index = obj.CornerPointGridData.Internal_Face.CellNeighbor2Index;
            CellNeighbor1Vec = obj.CornerPointGridData.Internal_Face.CellNeighbor1Vec;
            CellNeighbor2Vec = obj.CornerPointGridData.Internal_Face.CellNeighbor2Vec;
            Nvec = obj.CornerPointGridData.Internal_Face.Nvec;
            Trans_Half_1 = sum( Perm(CellNeighbor1Index,1) .* CellNeighbor1Vec .* Nvec , 2 ) ./ sum( CellNeighbor1Vec .* CellNeighbor1Vec , 2 );
            Trans_Half_2 = sum( Perm(CellNeighbor2Index,1) .* CellNeighbor2Vec .* Nvec , 2 ) ./ sum( CellNeighbor2Vec .* CellNeighbor2Vec , 2 );
            obj.Trans = Trans_Half_1 .* Trans_Half_2 ./ (Trans_Half_1 + Trans_Half_2);
        end
        function CorrectTransmissibilitiesForpEDFM(obj)
            % Virtual Call
        end
        function AddGridCoordinates(obj)
            obj.GridCoords = [obj.CornerPointGridData.Cell.NW_Top_Corner, obj.CornerPointGridData.Cell.NE_Top_Corner, ...
                              obj.CornerPointGridData.Cell.SW_Top_Corner, obj.CornerPointGridData.Cell.SE_Top_Corner, ...
                              obj.CornerPointGridData.Cell.NW_Bot_Corner, obj.CornerPointGridData.Cell.NE_Bot_Corner, ...
                              obj.CornerPointGridData.Cell.SW_Bot_Corner, obj.CornerPointGridData.Cell.SE_Bot_Corner];
        end
        function ComputeDepth(obj, alpha, Thickness)
            obj.Depth = obj.CornerPointGridData.Cell.Centroid(:,3);
        end
    end
end