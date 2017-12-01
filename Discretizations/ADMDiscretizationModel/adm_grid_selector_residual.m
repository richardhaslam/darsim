%  ADM Grid Selector: residual based
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 30 June 2017
%Last modified: 21 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid_selector_residual < adm_grid_selector
    properties
    end
    methods
        function obj = adm_grid_selector_residual(tol)
            obj@adm_grid_selector(tol);
        end
        function SelectGrid(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, Residual, maxLevel)
            % SELECT the ADM GRID for next time-step
            % Grid is chosen based on Residual(p^nu, S^nu), with nu = 1
            
            %% 1. Reset all cells to be active 
            n_media = length(FineGrid);
            Residuum = cell(n_media);
            Stop = 0;
            Ntot = sum(FineGrid(:).N);
            for m=1:n_media
                FineGrid(m).Active = ones(FineGrid(m).N, 1);
                if m==1
                    % for now wells are only in the reservoir
                    CoarseGrid(m,1).Active = obj.NoWellsCoarseCells;  
                else
                    CoarseGrid(m,1).Active = ones(CoarseGrid(m).N, 1);
                end
                Start = Stop + 1; 
                Stop = Start + FineGrid(m).N - 1;
                AbsResidual = [abs(Residual(Start:Stop)), abs(Residual(Ntot+Start:Ntot+Stop))];
                ScaledResidual = AbsResidual ./ max(AbsResidual);
                Residuum{m} = max(ScaledResidual(:, 1), ScaledResidual(:, 2));
            end
            figure(1)
            surf(reshape(Residuum{1}, 99, 99));
            drawnow
            
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
                        obj.SelectCoarseFine(FineGrid(m), CoarseGrid(m, 1), Residuum{m});
                    elseif l <= maxLevel(m)
                        obj.SelectCoarseFine(CoarseGrid(m, l-1), CoarseGrid(m, l), Residuum{m});
                    else
                        CoarseGrid(m, l).Active = zeros(CoarseGrid(m, l).N, 1);
                    end
                end
            end
            
            %% 3. Create ADM Grid
            obj.CreateADMGrid(ADMGrid, FineGrid, CoarseGrid, maxLevel);
        end
        function SelectCoarseFine(obj, FineGrid, CoarseGrid, Residual)
            %Given a Fine (level l-1) and a Coarse (level l) Grids chooses the cells that have to be active
            
            %% 2. Select Active Coarse Blocks
            Nc = CoarseGrid.N;
            for c = 1:Nc
                % fine-scale cells inside coarse block c
                indexes_fs = CoarseGrid.GrandChildren(c,:);
                
                % Max delta inside block c
                ResidualSum = max(Residual(indexes_fs)); 
   
                if CoarseGrid.Active(c) == 1 && abs(ResidualSum) > obj.tol
                   CoarseGrid.Active(c) = 0;
                end
            end
            
            %% 3. Set to inactive fine blocks (level l-1) belonging to active Coarse Blocks (level l)
            FineGrid.Active(CoarseGrid.Children(CoarseGrid.Active == 1,:)) = 0;
        end
    end
end