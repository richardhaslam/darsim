%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef RefCellsSelector < handle
    properties
        tol
        NSten = 2
        ActCells
        BCFluxMatrix
        ThsBC
        f
        Pc
        Mob2
        rho
    end
    methods
        function obj = RefCellsSelector(tol)
            obj.tol = tol;
        end
        function CopyCellsSelected(obj, CellsSelectedOld)
            % to create a vector of Cells Selected inside the sub-ref I
            % need to copy the class
            obj.ActCells = CellsSelectedOld.ActCells;
            obj.f = CellsSelectedOld.f;
            obj.Mob2 = CellsSelectedOld.Mob2;
            obj.Pc = CellsSelectedOld.Pc;
            obj.rho = CellsSelectedOld.rho;
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
            obj.f = Formulation.f;
            obj.Pc = ProductionSystem.Reservoir.State.Properties('Pc').Value;
            obj.Mob2 = Formulation.Mob(1:N, 2);
            obj.rho = ProductionSystem.Reservoir.State.Properties('rho_1').Value;
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
            
            MatrixAssembler.ActCells = obj.ActCells;
        end
     
        function ComputeActiveCells(obj, ProductionSystem, DiscretizationModel, Formulation, level)
            obj.ActCells = DiscretizationModel.FineGrid.Active;
            N = DiscretizationModel.ReservoirGrid.N;
            
            for i = level:  DiscretizationModel.maxLevel - 1
                obj.ActCells([DiscretizationModel.CoarseGrid(i).GrandChildren{DiscretizationModel.CoarseGrid(i).Active == 1,:}]) = 1;
            end
            obj.f = Formulation.f;
            obj.Pc = ProductionSystem.Reservoir.State.Properties('Pc').Value;
            obj.Mob2 = Formulation.Mob(1:N, 2);
            obj.rho = ProductionSystem.Reservoir.State.Properties('rho_1').Value;
            
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