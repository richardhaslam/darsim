%  ADM Grid Selector: delta prop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 30 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid_selector_delta < adm_grid_selector
    properties
    end
    methods
        function obj = adm_grid_selector_delta(tol, key)
            obj@adm_grid_selector(tol);
            obj.key = key;
        end
        function SelectGrid(obj, FineGrids, CoarseGrids, ADMGrid, ProductionSystem, Residual, maxLevel)
            % SELECT the ADM GRID for next time-step based on delta x
            %% 1. Reset all cells to be active and store Var{m}
            n_media = length(FineGrids);
            Var = cell(n_media, 1);
            for m=1:n_media
                FineGrids(m).Active(:) = 1;
                % Setting the region around the wells active as fine-scale resolution
                if m==1 % for now wells are only in the reservoir
                    CoarseGrids(m,1).Active = obj.NoWellsCoarseCells;
                    Var{m} = ProductionSystem.Reservoir.State.Properties(obj.key).Value;
                else 
                    CoarseGrids(m,1).Active(:) = 1;
                    Var{m} = ProductionSystem.FracturesNetwork.Fractures(m-1).State.Properties(obj.key).Value;
                end
            end
            
            %% 2. Select active cells
            for L=1:maxLevel(1)
                % 2.a choose possible active grids for level l
                if L>1
                    obj.DefinePossibleActive(CoarseGrids(:, L), CoarseGrids(:, L-1), L);
                end
                for m=1:n_media
                    % 2.b choose active cells of level l
                    if L==1
                        % coarse grid 1 to 0 (fine-scale)
                        obj.SelectCoarseFine(FineGrids(:,1), CoarseGrids(:,1), m, Var, false);
                    elseif L <= max(maxLevel)%maxLevel(m)
                        obj.SelectCoarseFine(CoarseGrids(:,L-1), CoarseGrids(:,L), m, Var, true);
                    else
                        CoarseGrids(m, L).Active = zeros(CoarseGrids(m,L).N, 1);
                    end
                end
            end
            %% 3. Create ADM Grid
            obj.CreateADMGrid(ADMGrid, FineGrids, CoarseGrids, maxLevel);
        end
        function SelectGridCoarse(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, Residual, maxLevel)
            % SELECT the ADM GRID for next time-step based on delta x
            %% 1. Reset all cells to be active and stor x{m}
            n_media = length(FineGrid);
            for m=1:n_media
                FineGrid(m).Active = ones(FineGrid(m).N, 1);
                if m==1
                    % for now wells are only in the reservoir
                    CoarseGrid(m,1).Active = obj.NoWellsCoarseCells;
                else
                    CoarseGrid(m,1).Active = ones(CoarseGrid(m).N, 1);
                end
            end
            
            %% 2. Select active cells
            for l=1:max(maxLevel)
                % 2.a choose possible active grids for level l
                if l>1
                    obj.DefinePossibleActive(CoarseGrid(:, l), CoarseGrid(:, l-1), l);
                end
                for m=1:n_media
                    % 2.b choose active cells of level l
                    if l==1
                        % coarse grid 1 to 0 (fine-scale)
                        obj.SelectCoarseFineC(FineGrid(m), CoarseGrid(m, 1));
                    elseif l <= max(maxLevel)%maxLevel(m)
                        obj.SelectCoarseFineC(CoarseGrid(m, l-1), CoarseGrid(m, l));
                    else
                        CoarseGrid(m, l).Active = zeros(CoarseGrid(m, l).N, 1);
                    end
                end
            end
            obj.CreateADMGrid(ADMGrid, FineGrid, CoarseGrid, maxLevel);
        end
        function SelectCoarseFine(obj, FineGrids, CoarseGrids, m, Var, flag)
            %Given a Fine and a Coarse Grid chooses the cells that have to be active
            
            %1. Select Active Coarse Blocks
            Nc = [CoarseGrids.N]';
            for c = 1:Nc(m)
                % Values of the variable "Var" for FineGrid cells belonging to CoarseGrid block c
                Var_children = Var{m}(CoarseGrids(m).Children{c, end});
                % Finding max and min values inside CoarseGrid block c
                VarMax = max(Var_children);
                VarMin = min(Var_children);
                if CoarseGrids(m).Active(c) == 1
                    % 1.1. Checking the neighbours (inside this medium itself)
                    Neighbour = CoarseGrids(m).Neighbours{c};
                    NofNeighbour = length(Neighbour);
                    i = 1;
                    while i <= NofNeighbour
                        % Finding max and min values of neighbour i
                        Var_children_Neighbor = Var{m}(CoarseGrids(m).Children{Neighbour(i), end});
                        VarNeighborMax = max(Var_children_Neighbor);
                        VarNeighborMin = min(Var_children_Neighbor);
                        if flag == false
                            if (abs(VarMax-VarNeighborMin) > obj.tol || abs(VarMin-VarNeighborMax) > obj.tol)
                                CoarseGrids(m).Active(c) = 0;
                                i = NofNeighbour + 1;
                            else
                                i = i+1;
                            end
                        else
                            if (abs(VarMax-VarNeighborMin) > 2*obj.tol || abs(VarMin-VarNeighborMax) > 2*obj.tol)
                                CoarseGrids(m).Active(c) = 0;
                                i = NofNeighbour + 1;
                            else
                                i = i+1;
                            end
                        end
                    end
                    
                    % 1.2. Checking the non-neighbours (inside other media) if any
                    if obj.isCoupled && length(CoarseGrids) > 1 % which means we have fractures
                        NonNeighbour = CoarseGrids(m).NonNeighbours{c};
                        NOfNonNeighbour = length(NonNeighbour);
                        i = 1;
                        while i <= NOfNonNeighbour
                            % Finding max and min values of non-neighbour i
                            Ind = cumsum(Nc) - NonNeighbour(i);
                            Ind(Ind<0) = nan;
                            [~, m_other] = min(Ind);
                            if m == m_other,  error("The non-neighboring connection cannot be from the same medium. Something is wrong.");  end
                            c_other = NonNeighbour(i) - sum(Nc(1:m_other-1));
                            Var_children_NonNeighbor = Var{m_other}(CoarseGrids(m_other).Children{c_other, end});
                            VarNonNeighborMax = max(Var_children_NonNeighbor);
                            VarNonNeighborMin = min(Var_children_NonNeighbor);
                            if flag == false
                                if (abs(VarMax-VarNonNeighborMin) > obj.tol || abs(VarMin-VarNonNeighborMax) > obj.tol)
                                    CoarseGrids(m).Active(c) = 0;
                                    i = NOfNonNeighbour + 1;
                                else
                                    i = i+1;
                                end
                            else
                                if (abs(VarMax-VarNonNeighborMin) > 2*obj.tol || abs(VarMin-VarNonNeighborMax) > 2*obj.tol)
                                    CoarseGrids(m).Active(c) = 0;
                                    i = NOfNonNeighbour + 1;
                                else
                                    i = i+1;
                                end
                            end
                        end
                    end
                end
            end
            
            %2. Do not coarsen neighbours of cells that are at finescale resolution
            DummyActive = CoarseGrids(m).Active;
            for i = 1:Nc(m)
                if (CoarseGrids(m).Active(i) == 0)
                    DummyActive(CoarseGrids(m).Neighbours{i}) = 0;
                end
            end
            CoarseGrids(m).Active = DummyActive.*CoarseGrids(m).Active;
            
            %3. Set to inactive the fine block belonging to Active Coarse Blocks
            FineGrids(m).Active([CoarseGrids(m).Children{CoarseGrids(m).Active == 1,1}]) = 0;
        end
        function SelectCoarseFineC(obj, FineGrid, CoarseGrid)
            FineGrid.Active([CoarseGrid.Children{CoarseGrid.Active == 1,:}]) = 0;
        end
    end
end