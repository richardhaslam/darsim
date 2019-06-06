%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef RefCellsSelector < handle
    properties
        tol
        NSten = 2
        ActFluxes
        ActCells
        Dirichlet
        ActCellsMesh
        ViscousMatrixValue
        f
    end
    methods
        function obj = RefCellsSelector(tol)
            obj.tol = tol;
        end
        function CopyCellsSelected(obj, CellsSelectedOld)
            % to create a vector of Cells Selected inside the sub-ref I
            % need to copy the class
            obj.ActFluxes = CellsSelectedOld.ActFluxes;
            obj.ActCells = CellsSelectedOld.ActCells;
            obj.ViscousMatrixValue = CellsSelectedOld.ViscousMatrixValue;
            obj.f = CellsSelectedOld.f;
        end
        function SelectRefCells(obj, ProductionSystem, Grid, Formulation)
            
            Utot = Formulation.Utot;
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            % We choose to resolve with a fine dt all cells which were
            % flooded during the coarse dt.
            % For now it only works for the reservoir and not for the
            % fractures
            Snew = ProductionSystem.Reservoir.State.Properties('S_1').Value;
            Sold = ProductionSystem.Reservoir.State_old.Properties('S_1').Value;
            dS = Snew - Sold;
            % vector of active cells
            obj.ActCells = abs(dS) >= obj.tol;
            % obj.ActCells = (Snew < 0.65 & Snew > 0.105);
            
            % to add also the next cell according to the upwind velocity
            
            ActCellsM = reshape(obj.ActCells, Nx, Ny, Nz);
            
            % if you want to increase the set just apply a for cycle
            for i = 1: obj.NSten
                ActCellsM(2:Nx,:,:) = ActCellsM(2:Nx,:,:) + max(Utot.x(2:Nx,:,:),0) .* ActCellsM(1:Nx-1,:,:);
                ActCellsM(1:Nx-1,:,:) = ActCellsM(1:Nx-1,:,:) + min(Utot.x(1:Nx-1,:,:),0) .*  ActCellsM(2:Nx,:,:);
                
                ActCellsM(:,2:Ny,:) = ActCellsM(:,2:Ny,:) + max(Utot.y(:,2:Ny,:),0) .* ActCellsM(:,1:Ny-1,:);
                ActCellsM(:,1:Ny-1,:) = ActCellsM(:,1:Ny-1,:) + min(Utot.y(:,1:Ny-1,:),0) .* ActCellsM(:,2:Ny,:);
                
                ActCellsM(:,:,2:Nz) = ActCellsM(:,:,2:Nz) + max(Utot.z(:,:,2:Nz),0) .* ActCellsM(:,:,1:Nz-1);
                ActCellsM(:,:,1:Nz-1) = ActCellsM(:,:,1:Nz-1) + min(Utot.z(:,:,1:Nz-1),0) .* ActCellsM(:,:,2:Nz);
            end
            obj.ActCells = reshape(ActCellsM >0, N, 1);
        end
        function SetActiveInterfaces(obj, MatrixAssembler, Grid)
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            
            % From the ActCells mask we define the Active Fluxes.
            MatrixAssembler.ActInterfaces.x = zeros(Nx+1,Ny,Nz);
            MatrixAssembler.ActInterfaces.y = zeros(Nx,Ny+1,Nz);
            MatrixAssembler.ActInterfaces.z = zeros(Nx,Ny,Nz+1);
            
            ActCellsM = reshape(obj.ActCells, Nx, Ny, Nz);
            
            MatrixAssembler.ActInterfaces.x(1,:,:)    = ActCellsM(1,:,:);
            MatrixAssembler.ActInterfaces.x(Nx+1,:,:) = ActCellsM(Nx,:,:);
            MatrixAssembler.ActInterfaces.x(2:Nx,:,:) = (ActCellsM(1:Nx-1,:,:) + ActCellsM(2:Nx,:,:)) == 2;
            
            MatrixAssembler.ActInterfaces.y(:,1,:)    = ActCellsM(:,1,:);
            MatrixAssembler.ActInterfaces.y(:,Ny+1,:) = ActCellsM(:,Ny,:);
            MatrixAssembler.ActInterfaces.y(:,2:Ny,:) = (ActCellsM(:,1:Ny-1,:) + ActCellsM(:,2:Ny,:)) == 2;
            
            MatrixAssembler.ActInterfaces.z(:,:,1)    = ActCellsM(:,:,1);
            MatrixAssembler.ActInterfaces.z(:,:,Nz+1) = ActCellsM(:,:,Nz);
            MatrixAssembler.ActInterfaces.z(:,:,2:Nz) = (ActCellsM(:,:,1:Nz-1) + ActCellsM(:,:,2:Nz)) == 2;
        end
        function ActFluxes = SelectRefFluxes(obj, Grid)
            
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            
            % From the ActCells mask we define the Active Fluxes.
            ActFluxes.x = zeros(Nx+1,Ny,Nz);
            ActFluxes.y = zeros(Nx,Ny+1,Nz);
            ActFluxes.z = zeros(Nx,Ny,Nz+1);
            
            ActCellsM = reshape(obj.ActCells, Nx, Ny, Nz);
            
            ActFluxes.x(1,:,:)    = ActCellsM(1,:,:);
            ActFluxes.x(Nx+1,:,:) = ActCellsM(Nx,:,:);
            ActFluxes.x(2:Nx,:,:) = (ActCellsM(1:Nx-1,:,:) + ActCellsM(2:Nx,:,:)) == 2;
            
            ActFluxes.y(:,1,:)    = ActCellsM(:,1,:);
            ActFluxes.y(:,Ny+1,:) = ActCellsM(:,Ny,:);
            ActFluxes.y(:,2:Ny,:) = (ActCellsM(:,1:Ny-1,:) + ActCellsM(:,2:Ny,:)) == 2;
            
            ActFluxes.z(:,:,1)    = ActCellsM(:,:,1);
            ActFluxes.z(:,:,Nz+1) = ActCellsM(:,:,Nz);
            ActFluxes.z(:,:,2:Nz) = (ActCellsM(:,:,1:Nz-1) + ActCellsM(:,:,2:Nz)) == 2;
            obj.ActFluxes = ActFluxes;
        end
        function ComputeBoundaryValues(obj, DiscretizationModel, Formulation)
            % store the past values of the fractional flow
            obj.f = Formulation.f;
            
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N  = DiscretizationModel.ReservoirGrid.N;
            
            BCFluxes.x = zeros(Nx+1,Ny,Nz);
            BCFluxes.y = zeros(Nx,Ny+1,Nz);
            BCFluxes.z = zeros(Nx,Ny,Nz+1);
            ActCellsM = reshape(obj.ActCells, Nx, Ny, Nz);
            
            % Compute the fluxes values
            Utot = Formulation.Utot;
            
            % right to left and top to bottom (negative x, y, z)
            Xneg = min(Utot.x, 0);
            Yneg = min(Utot.y, 0);
            Zneg = min(Utot.z, 0);
            
            x1 = reshape(Xneg(1:Nx,:,:),N,1);
            y1 = reshape(Yneg(:,1:Ny,:),N,1);
            z1 = reshape(Zneg(:,:,1:Nz),N,1);
            
            % left to right and bottom to top (positive x, y, z)
            Xpos = max(Utot.x, 0);
            Ypos = max(Utot.y, 0);
            Zpos = max(Utot.z, 0);
            
            x2 = reshape(Xpos(2:Nx+1,:,:), N, 1);
            y2 = reshape(Ypos(:,2:Ny+1,:), N, 1);
            z2 = reshape(Zpos(:,:,2:Nz+1), N, 1);
            
            
            %fluxes at the boundaty between an accepted and a refuced cell
            BCFluxes.x(2:Nx,:,:) = (ActCellsM(1:Nx-1,:,:) + ActCellsM(2:Nx,:,:) == 1);
            BCFluxes.y(:,2:Ny,:) = (ActCellsM(:,1:Ny-1,:) + ActCellsM(:,2:Ny,:) == 1);
            BCFluxes.z(:,:,2:Nz) = (ActCellsM(:,:,1:Nz-1) + ActCellsM(:,:,2:Nz) == 1);
            
            
            %select only the bc fluxes
            x1A = reshape(BCFluxes.x(1:Nx,:,:),N,1);
            y1A = reshape(BCFluxes.y(:,1:Ny,:),N,1);
            z1A = reshape(BCFluxes.z(:,:,1:Nz),N,1);
            
            x1 = x1 .* x1A;
            y1 = y1 .* y1A;
            z1 = z1 .* z1A;
            
            %select only the bc fluxes
            x2A = reshape(BCFluxes.x(2:Nx+1,:,:),N,1);
            y2A = reshape(BCFluxes.y(:,2:Ny+1,:),N,1);
            z2A = reshape(BCFluxes.z(:,:,2:Nz+1),N,1);
            
            x2 = x2 .* x2A;
            y2 = y2 .* y2A;
            z2 = z2 .* z2A;
            
            DiagVecs = [z2, y2, x2, +x1-x2+y1-y2+z1-z2, -x1, -y1, -z1]; % diagonal vectors
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            
            obj.ViscousMatrixValue = spdiags(DiagVecs, DiagIndx, N, N);
        end
        function ComputeBoundaryValuesSubRef(obj, DiscretizationModel, Formulation, CellsSelected_old)
            % merge the new accepted fractional flow with those occepted at
            % the previous sub-ref
            
            obj.f = Formulation.f .* (CellsSelected_old.ActCells) + ...
                CellsSelected_old.f .* (1 - CellsSelected_old.ActCells);
            
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N  = DiscretizationModel.ReservoirGrid.N;
            % merge the new fluxes interfaces with those obtained at the
            % previous sub-ref.
            
            BCFluxes.x = zeros(Nx+1,Ny,Nz);
            BCFluxes.y = zeros(Nx,Ny+1,Nz);
            BCFluxes.z = zeros(Nx,Ny,Nz+1);
            ActCellsM = reshape(obj.ActCells, Nx, Ny, Nz);
            % Compute the fluxes values
            Utot = Formulation.Utot;
            
            % right to left and top to bottom (negative x, y, z)
            Xneg = min(Utot.x, 0);
            Yneg = min(Utot.y, 0);
            Zneg = min(Utot.z, 0);
            
            x1 = reshape(Xneg(1:Nx,:,:),N,1);
            y1 = reshape(Yneg(:,1:Ny,:),N,1);
            z1 = reshape(Zneg(:,:,1:Nz),N,1);
            
            % left to right and bottom to top (positive x, y, z)
            Xpos = max(Utot.x, 0);
            Ypos = max(Utot.y, 0);
            Zpos = max(Utot.z, 0);
            
            x2 = reshape(Xpos(2:Nx+1,:,:), N, 1);
            y2 = reshape(Ypos(:,2:Ny+1,:), N, 1);
            z2 = reshape(Zpos(:,:,2:Nz+1), N, 1);
            
            %fluxes at the boundary between an accepted and a rejected cell
            BCFluxes.x(2:Nx,:,:) = (ActCellsM(1:Nx-1,:,:) + ActCellsM(2:Nx,:,:) == 1);
            BCFluxes.y(:,2:Ny,:) = (ActCellsM(:,1:Ny-1,:) + ActCellsM(:,2:Ny,:) == 1);
            BCFluxes.z(:,:,2:Nz) = (ActCellsM(:,:,1:Nz-1) + ActCellsM(:,:,2:Nz) == 1);
            
            
            %select only the bc fluxes
            x1A = reshape(BCFluxes.x(1:Nx,:,:),N,1);
            y1A = reshape(BCFluxes.y(:,1:Ny,:),N,1);
            z1A = reshape(BCFluxes.z(:,:,1:Nz),N,1);
            
            x1 = x1 .* x1A;
            y1 = y1 .* y1A;
            z1 = z1 .* z1A;
            
            %select only the bc fluxes
            x2A = reshape(BCFluxes.x(2:Nx+1,:,:),N,1);
            y2A = reshape(BCFluxes.y(:,2:Ny+1,:),N,1);
            z2A = reshape(BCFluxes.z(:,:,2:Nz+1),N,1);
            
            x2 = x2 .* x2A;
            y2 = y2 .* y2A;
            z2 = z2 .* z2A;
            
            DiagVecs = [z2, y2, x2, +x1-x2+y1-y2+z1-z2, -x1, -y1, -z1]; % diagonal vectors
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            
            obj.ViscousMatrixValue = spdiags(DiagVecs, DiagIndx, N, N);
        end    
        function ComputeActiveCells(obj, DiscretizationModel, level)
            obj.ActCells = DiscretizationModel.FineGrid.Active;
            for i = level:  DiscretizationModel.maxLevel - 1
                obj.ActCells([DiscretizationModel.CoarseGrid(i).GrandChildren{DiscretizationModel.CoarseGrid(i).Active == 1,:}]) = 1;
            end
        end
        function numb = NumberOfActiveCells(obj, DiscretizationModel, level)
            if level == DiscretizationModel.maxLevel
                numb = sum(obj.ActCells);
            else
                numb = sum(DiscretizationModel.FineGrid.Active .* obj.ActCells);
                for i = level:  DiscretizationModel.maxLevel - 1
                    CellCoarse = zeros(size(obj.ActCells));
                    % I need just the corrispondence between one cell of
                    % the fine that belongs to the coarse to count the
                    % number of Active in time coarse cells
                    CellCoarse ([DiscretizationModel.CoarseGrid(i).GrandChildren{DiscretizationModel.CoarseGrid(i).Active == 1,1}]) = 1;
                    numb = numb + sum(CellCoarse .* obj.ActCells);
                end
            end
        end
        function ActCellCheckError(obj, ProductionSystem, Grid, Formulation)
            Utot = Formulation.Utot;
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            dS = ProductionSystem.Reservoir.State.Properties('S_1').Value - ...
                 ProductionSystem.Reservoir.State_old.Properties('S_1').Value;
            % vector of active cells
            PossibleCells = obj.ActCells;
            obj.ActCells = (abs(dS) >= obj.tol) .* PossibleCells ;
            ActCellsM = reshape(obj.ActCells, Nx, Ny, Nz);
            for i = 1: obj.NSten
                ActCellsM(2:Nx,:,:) = ActCellsM(2:Nx,:,:) + max(Utot.x(2:Nx,:,:),0) .* ActCellsM(1:Nx-1,:,:);
                ActCellsM(1:Nx-1,:,:) = ActCellsM(1:Nx-1,:,:) + min(Utot.x(1:Nx-1,:,:),0) .*  ActCellsM(2:Nx,:,:);
                
                ActCellsM(:,2:Ny,:) = ActCellsM(:,2:Ny,:) + max(Utot.y(:,2:Ny,:),0) .* ActCellsM(:,1:Ny-1,:);
                ActCellsM(:,1:Ny-1,:) = ActCellsM(:,1:Ny-1,:) + min(Utot.y(:,1:Ny-1,:),0) .* ActCellsM(:,2:Ny,:);
                
                ActCellsM(:,:,2:Nz) = ActCellsM(:,:,2:Nz) + max(Utot.z(:,:,2:Nz),0) .* ActCellsM(:,:,1:Nz-1);
                ActCellsM(:,:,1:Nz-1) = ActCellsM(:,:,1:Nz-1) + min(Utot.z(:,:,1:Nz-1),0) .* ActCellsM(:,:,2:Nz);
            end
            obj.ActCells = PossibleCells  .* reshape(ActCellsM >0, N, 1);
        end
        function ActiveCellsOnLevelAndNeighbours(obj, DiscretizationModel, level)
            %level = DiscretizationModel.maxLevel - level + 1;
            obj.ActCells= zeros(size(DiscretizationModel.FineGrid.Active));
            CoarseGrid = DiscretizationModel.CoarseGrid(level);
            ActCellCoarse = CoarseGrid.Active;
            Nc = CoarseGrid.N;
            for c = 1:Nc
                if CoarseGrid.Active(c) == 1
                    n = CoarseGrid.Neighbours(c).indexes;
                    Nn = length(n);
                    i = 1;
                    while i <= Nn
                        ActCellCoarse(n(i)) = 1;
                        i = i + 1;
                    end
                end
            end
            obj.ActCells(DiscretizationModel.CoarseGrid(level).GrandChildren(ActCellCoarse == 1,:)) = 1;
            ActLevel = zeros(size(DiscretizationModel.FineGrid.Active));
            ActLevel(DiscretizationModel.CoarseGrid(level).GrandChildren(DiscretizationModel.CoarseGrid(level).Active == 1,:)) = 1;
            if level > 1
                ActLevel(DiscretizationModel.CoarseGrid(level-1).GrandChildren(DiscretizationModel.CoarseGrid(level-1).Active == 1,:)) = 1;
            else
                ActLevel(DiscretizationModel.FineGrid.Active(:)==1) = 1;
            end
                obj.ActCells = obj.ActCells .* ActLevel;
        end
        function EffectiveActiveCellsDD(obj, DiscretizationModel, level)
            obj.ActCells= zeros(size(DiscretizationModel.FineGrid.Active));
            CoarseGrid = DiscretizationModel.CoarseGrid(level);
            ActCellCoarse = CoarseGrid.Active;
            obj.ActCells(DiscretizationModel.CoarseGrid(level).GrandChildren(ActCellCoarse == 1,:)) = 1;
        end
        function ResetInitialState(obj, State, State_old)
            Names = State.Properties.keys;
            N_prop = double(State.Properties.Count);
            for i=1:N_prop
                temp = State.Properties(Names{i});
                temp.Value(obj.ActCells) = State_old.Properties(Names{i}).Value(obj.ActCells);
            end
        end
    end
end