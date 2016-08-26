%  ADM operators handler base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 16 August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_handler_MS < operators_handler
    properties
        BFUpdater
    end
    methods
        function obj = operators_handler_MS(n)
            obj@operators_handler(n)
        end
        function BuildStaticOperators(obj, CoarseGrid, FineGrid, maxLevel, K, s, FluidModel)
            obj.BFUpdater.ConstructPressureSystem(FineGrid, K, s, FluidModel);
            %Build static restriction operator (FV)
            obj.R{1} = obj.MsRestriction(FineGrid, CoarseGrid(1), 1);
            % Build Prolongation operator
            obj.Pp{1} = obj.BFUpdater.MsProlongation(FineGrid, CoarseGrid(1), CoarseGrid(1).CoarseFactor);
            %Build first coarse system (with MsFE)
            obj.BFUpdater.A = obj.Pp{1}' * obj.BFUpdater.A * obj.Pp{1};
            for x = 2:maxLevel
                % Build static restriction operator (FV)
                obj.R{x} = obj.MsRestriction(CoarseGrid(x-1), CoarseGrid(x), x);
                % Build Prolongation operator
                obj.Pp{x} = obj.BFUpdater.MsProlongation(CoarseGrid(x-1), CoarseGrid(x), CoarseGrid(x).CoarseFactor./CoarseGrid(x-1).CoarseFactor);
                %Build coarse system (with MsFE)
                obj.BFUpdater.A = obj.Pp{x}' * obj.BFUpdater.A * obj.Pp{x};
            end
        end
        function [R, Pp, Ps] = BuildADMOperators(FineGrid, CoarseGrid, ADMGrid)
            %Construct DLGR R and P
            maxLevel = ADMGrid.level(end);
            %Other Levels
            for i=maxLevel:-1:2
                [Nf, Nx] = NewNumberOfCells(DLGRGrid, i);
                [R(i).matrix, Pp(i).matrix, Ps(i).matrix, ADMGrid] = BuildRandP(CoarseGrid(i-1), CoarseGrid(i), ADMGrid, i, sum(ADMGrid.N), Nf, Nx);
            end
            %First level
            [R(1).matrix, Pp(1).matrix, Ps(1).matrix] = BuildRandP(FineGrid, CoarseGrid(1), ADMGrid, 1, sum(ADMGrid.N), FineGrid.Nx*FineGrid.Ny, ADMGrid.N(1));
        end
        
        function [Nf, Nx] = NewNumberOfCells(ADMGrid, x)
            %For now I use 9 that is CF between two levels
            Nx = sum(ADMGrid.N) - ADMGrid.N(x+1);
            Nf = Nx + ADMGrid.N(x+1) * 9;
        end
    end
end