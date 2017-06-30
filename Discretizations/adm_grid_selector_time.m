%  ADM Grid Selector: time criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 30 June 2017
%Last modified: 30 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid_selector_time < adm_grid_selector
    properties
    end
    methods
        function obj = adm_grid_selector_time(tol, key)
            obj@adm_grid_selector(tol, key);
        end
        function SelectGrid(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, maxLevel)
            % SELECT the ADM GRID for next time-step
            FineGrid.Active = ones(FineGrid.N, 1);
           
            obj.CreateADMGrid(ADMGrid, FineGrid, CoarseGrid);
        end
        function SelectCoarseFine(obj, FineGrid, CoarseGrid, S)
            %Given a Fine and a Coarse Grids chooses the cells that have to be active
            %1. Select Active Coarse Blocks
            Nc = CoarseGrid.N;
            for c = 1:Nc
                I = CoarseGrid.I(c,2);
                J = CoarseGrid.J(c,2);
                K = CoarseGrid.K(c,2);
                Imin = I - floor((CoarseGrid.CoarseFactor(1) - 1)/2);
                Imax = I + ceil((CoarseGrid.CoarseFactor(1) - 1)/2);
                Jmin = J - floor((CoarseGrid.CoarseFactor(2) - 1)/2);
                Jmax = J + ceil((CoarseGrid.CoarseFactor(2) - 1)/2);
                Kmin = K - floor((CoarseGrid.CoarseFactor(3) - 1)/2);
                Kmax = K + ceil((CoarseGrid.CoarseFactor(3) - 1)/2);
                
                % Max e Min saturation
                Smax = max(max(max(S(Imin:Imax, Jmin:Jmax, Kmin:Kmax))));
                Smin = min(min(min(S(Imin:Imax, Jmin:Jmax, Kmin:Kmax))));
                if CoarseGrid.Active(c) == 1
                    n = CoarseGrid.Neighbours(c).indexes;
                    Nn = length(n);
                    i = 1;
                    while i <= Nn
                        if (abs(Smax-S(CoarseGrid.I(n(i), 2), CoarseGrid.J(n(i), 2), CoarseGrid.K(n(i), 2)))...
                                > obj.tol || abs(Smin-S(CoarseGrid.I(n(i),2),CoarseGrid.J(n(i),2), CoarseGrid.K(n(i), 2))) > obj.tol)
                            CoarseGrid.Active(c) = 0;
                            %CoarseGrid.Active(i) = 0;
                            i = Nn + 1;
                        else
                            i = i+1;
                        end
                    end
                end
            end
            
            %2. Do not coarsen neighbours of cells that are fine
            DummyActive = CoarseGrid.Active;
            for i = 1:Nc
                if (CoarseGrid.Active(i) == 0)
                    DummyActive(CoarseGrid.Neighbours(i).indexes) = 0;
                end
            end
            CoarseGrid.Active = DummyActive.*CoarseGrid.Active;
            
            %3. Set to inactive fine block belonging to Active Coarse Blocks
            %Cindeces = find();
            FineGrid.Active(CoarseGrid.Children(CoarseGrid.Active == 1,:)) = 0;
        end
    end
end