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
        pEDFM_alpha_Trans
        HeatTrans
        Nx
        Ny
        Nz
        Volume
        ConnectivityMatrix
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
            obj.pEDFM_alpha_Trans = zeros(ReservoirProperties.CornerPointGridData.N_InternalFaces,1);
            obj.HeatTrans = zeros(ReservoirProperties.CornerPointGridData.N_InternalFaces,1);
        end
        function Initialize(obj, Reservoir)
            obj.Volume = obj.CornerPointGridData.Cell.Volume;
            obj.Neighbours = obj.CornerPointGridData.Cell.Index_Neighbors;
            obj.ComputeRockTransmissibilities(Reservoir.K);
            if ~isempty(Reservoir.K_Cond_eff)
                obj.ComputeRockHeatConductivities(Reservoir.K_Cond_eff)
            end
            obj.ConstructConnectivityMatrix();
        end
        function ComputeRockTransmissibilities(obj, Permeability)
            CellNeighbor1Index = obj.CornerPointGridData.Internal_Face.CellNeighbor1Index;
            CellNeighbor2Index = obj.CornerPointGridData.Internal_Face.CellNeighbor2Index;
            CellNeighbor1Vec = obj.CornerPointGridData.Internal_Face.CellNeighbor1Vec;
            CellNeighbor2Vec = obj.CornerPointGridData.Internal_Face.CellNeighbor2Vec;
            Nvec = obj.CornerPointGridData.Internal_Face.Nvec;
            Trans_Half_1 = sum( Permeability(CellNeighbor1Index,1) .* CellNeighbor1Vec .* Nvec , 2 ) ./ sum( CellNeighbor1Vec .* CellNeighbor1Vec , 2 );
            Trans_Half_2 = sum( Permeability(CellNeighbor2Index,1) .* CellNeighbor2Vec .* Nvec , 2 ) ./ sum( CellNeighbor2Vec .* CellNeighbor2Vec , 2 );
            Trans_Half_1 = abs(Trans_Half_1);
            Trans_Half_2 = abs(Trans_Half_2);
            obj.Trans = Trans_Half_1 .* Trans_Half_2 ./ (Trans_Half_1 + Trans_Half_2);
        end
        function AddpEDFMCorrections(obj,pEDFM_alpha_Trans)
            obj.pEDFM_alpha_Trans = pEDFM_alpha_Trans;
        end
        function CorrectTransmissibilitiesForpEDFM(obj)
            obj.Trans = obj.Trans .* ( 1 - obj.pEDFM_alpha_Trans );
        end
        function ComputeRockHeatConductivities(obj, K_Cond)
            CellNeighbor1Index = obj.CornerPointGridData.Internal_Face.CellNeighbor1Index;
            CellNeighbor2Index = obj.CornerPointGridData.Internal_Face.CellNeighbor2Index;
            CellNeighbor1Vec = obj.CornerPointGridData.Internal_Face.CellNeighbor1Vec;
            CellNeighbor2Vec = obj.CornerPointGridData.Internal_Face.CellNeighbor2Vec;
            Nvec = obj.CornerPointGridData.Internal_Face.Nvec;
            Trans_Half_1 = sum( K_Cond(CellNeighbor1Index,1) .* CellNeighbor1Vec .* Nvec , 2 ) ./ sum( CellNeighbor1Vec .* CellNeighbor1Vec , 2 );
            Trans_Half_2 = sum( K_Cond(CellNeighbor2Index,1) .* CellNeighbor2Vec .* Nvec , 2 ) ./ sum( CellNeighbor2Vec .* CellNeighbor2Vec , 2 );
            Trans_Half_1 = abs(Trans_Half_1);
            Trans_Half_2 = abs(Trans_Half_2);
            obj.HeatTrans = Trans_Half_1 .* Trans_Half_2 ./ (Trans_Half_1 + Trans_Half_2);
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
        function ConstructConnectivityMatrix(obj)
            % Creating a sparse matrix with number of rows being number of
            % the grid cells and number of columns being the number of
            % interfaces
            nc = obj.N;
            nf = length(obj.Trans);
            C = [ obj.CornerPointGridData.Internal_Face.CellNeighbor1Index , obj.CornerPointGridData.Internal_Face.CellNeighbor2Index ];
            obj.ConnectivityMatrix = sparse([(1:nf)'; (1:nf)'], C, ones(nf,1)*[-1 1], nf, nc)';
        end
    end
end