%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef RefCellsSelector < handle
    properties
        tol = 5e-2
        NSten = 2
        ActFluxes
        ActCells
        ActCellsMesh
        BCEnforcer = LTS_bc_enforcer_seq();
    end
    methods
        function obj = RefCellsSelector(coupling)
            switch(coupling)
                case('FIM')
                    obj.BCEnforcer = LTS_bc_enforcer_fim();
                case('Sequential')
                    obj.BCEnforcer = LTS_bc_enforcer_seq();
            end
        end
        function CopyCellsSelected(obj, CellsSelectedOld)
            % to create a vector of Cells Selected inside the sub-ref I
            % need to copy the class
            obj.ActFluxes = CellsSelectedOld.ActFluxes;
            obj.ActCells = CellsSelectedOld.ActCells;
            obj.BCEnforcer.ViscousMatrixValue = CellsSelectedOld.BCEnforcer.ViscousMatrixValue;
            obj.BCEnforcer.f = CellsSelectedOld.BCEnforcer.f;
        end
        function SelectRefCells(obj, ProductionSystem, Grid, Formulation)
            
            Utot = Formulation.Utot;
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            % We choose to resolve with a fine dt all cells which were
            % flooded during the coarse dt.
            %For now it only works for the reservoir and not for the
            %fractures
            dS = ProductionSystem.Reservoir.State.Properties('S_1').Value - ...
                ProductionSystem.Reservoir.State_old.Properties('S_1').Value;
            % vector of active cells
            obj.ActCells = abs(dS) >= obj.tol;
            
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
        function SetActiveInterfaces(obj, Formulation, Grid)
            Formulation.MatrixAssembler.ActiveInterfaces = ...
                obj.SelectRefFluxes(Grid);
            % Set velocities of inactive interfaces to be equal to zero.
            for ph = 1:2
                Formulation.U{ph, 1}.x = Formulation.U{ph, 1}.x .* obj.ActFluxes.x;
                Formulation.U{ph, 1}.y = Formulation.U{ph, 1}.y .* obj.ActFluxes.y;
                Formulation.U{ph, 1}.z = Formulation.U{ph, 1}.z .* obj.ActFluxes.z;
            end
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
            obj.BCEnforcer.ComputeBoundaryValues(DiscretizationModel, Formulation, obj.ActCells)
        end
        function ComputeBoundaryValuesSubRef(obj, DiscretizationModel, Formulation, CellsSelected_old)
            obj.BCEnforcer.ComputeBoundaryValuesSubRef(DiscretizationModel, Formulation, CellsSelected_old, obj.ActCells)
        end
        function ComputeActiveCells(obj, DiscretizationModel, level)
            obj.ActCells = DiscretizationModel.FineGrid.Active;
            for i = level:  DiscretizationModel.maxLevel - 1
                obj.ActCells(DiscretizationModel.CoarseGrid(i).GrandChildren(DiscretizationModel.CoarseGrid(i).Active == 1,:)) = 1;
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
                    CellCoarse (DiscretizationModel.CoarseGrid(i).GrandChildren(DiscretizationModel.CoarseGrid(i).Active == 1,1)) = 1;
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
    end
end