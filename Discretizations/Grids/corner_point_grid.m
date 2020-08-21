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
        dx
        dy
        dz
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
            obj.Volume = obj.CornerPointGridData.Cells.Volume;
            obj.Neighbours = obj.CornerPointGridData.Cells.Neighbors;
            obj.AddGridCellSize;
            obj.ComputeRockTransmissibilities(Reservoir.K);
            if ~isempty(Reservoir.K_Cond_eff)
                obj.ComputeRockHeatConductivities(Reservoir.K_Cond_eff)
            end
            obj.ConstructConnectivityMatrix();
        end
        function ComputeRockTransmissibilities(obj, Permeability)
            CellNeighbor1Index = obj.CornerPointGridData.Internal_Faces.CellNeighbor1Index;
            CellNeighbor2Index = obj.CornerPointGridData.Internal_Faces.CellNeighbor2Index;
            CellNeighbor1Vec = obj.CornerPointGridData.Internal_Faces.CellNeighbor1Vec;
            CellNeighbor2Vec = obj.CornerPointGridData.Internal_Faces.CellNeighbor2Vec;
            Nvec = obj.CornerPointGridData.Internal_Faces.Nvec;
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
            CellNeighbor1Index = obj.CornerPointGridData.Internal_Faces.CellNeighbor1Index;
            CellNeighbor2Index = obj.CornerPointGridData.Internal_Faces.CellNeighbor2Index;
            CellNeighbor1Vec = obj.CornerPointGridData.Internal_Faces.CellNeighbor1Vec;
            CellNeighbor2Vec = obj.CornerPointGridData.Internal_Faces.CellNeighbor2Vec;
            Nvec = obj.CornerPointGridData.Internal_Faces.Nvec;
            Trans_Half_1 = sum( K_Cond(CellNeighbor1Index,1) .* CellNeighbor1Vec .* Nvec , 2 ) ./ sum( CellNeighbor1Vec .* CellNeighbor1Vec , 2 );
            Trans_Half_2 = sum( K_Cond(CellNeighbor2Index,1) .* CellNeighbor2Vec .* Nvec , 2 ) ./ sum( CellNeighbor2Vec .* CellNeighbor2Vec , 2 );
            Trans_Half_1 = abs(Trans_Half_1);
            Trans_Half_2 = abs(Trans_Half_2);
            obj.HeatTrans = Trans_Half_1 .* Trans_Half_2 ./ (Trans_Half_1 + Trans_Half_2);
        end
        function AddGridCoordinates(obj)
            obj.GridCoords = obj.CornerPointGridData.Cells.Vertices;
        end
        function AddGridCellSize(obj)
            if sum(strcmp(fieldnames(obj.CornerPointGridData.Cells), 'dx'))
                obj.dx = obj.CornerPointGridData.Cells.dx;
                obj.dy = obj.CornerPointGridData.Cells.dy;
                obj.dz = obj.CornerPointGridData.Cells.dz;
            else
                obj.dx = zeros(obj.CornerPointGridData.N_ActiveCells,1);
                obj.dy = zeros(obj.CornerPointGridData.N_ActiveCells,1);
                obj.dz = zeros(obj.CornerPointGridData.N_ActiveCells,1);
                for i = 1 : obj.CornerPointGridData.N_ActiveCells
                    obj.dx(i) = max( obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(i,:),1) ) - ...
                                min( obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(i,:),1) );
                    obj.dy(i) = max( obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(i,:),2) ) - ...
                                min( obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(i,:),2) );
                    obj.dz(i) = max( obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(i,:),3) ) - ...
                                min( obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(i,:),3) );
                end
            end
        end
        function ComputeDepth(obj, alpha, Thickness)
            obj.Depth = max(obj.CornerPointGridData.Cells.Centroid(:,3)) - obj.CornerPointGridData.Cells.Centroid(:,3);
        end
        function ConstructConnectivityMatrix(obj)
            % Creating a sparse matrix with number of rows being number of
            % the grid cells and number of columns being the number of
            % interfaces
            nc = obj.N;
            nf = length(obj.Trans);
            C = [ obj.CornerPointGridData.Internal_Faces.CellNeighbor1Index , obj.CornerPointGridData.Internal_Faces.CellNeighbor2Index ];
            obj.ConnectivityMatrix = sparse([(1:nf)'; (1:nf)'], C, ones(nf,1)*[-1 1], nf, nc)';
        end
    end
end