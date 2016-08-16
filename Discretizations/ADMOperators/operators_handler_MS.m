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
        function BuildStaticOperators(obj, CoarseGrid, FineGrid, maxLevel, K)
            
            obj.BFUpdater.ConstructPressureSystem();
            %Build MS operators
            [CoarseGrid(1).MsR, CoarseGrid(1).MsP, CoarseGrid(1).C] = MSFVOperators(FineGrid, CoarseGrid(1), Ap, 1);
            %Build first coarse system (with MsFV)
            CoarseGrid(1).A_c = CoarseGrid(1).MsR*Ap*CoarseGrid(1).MsP;
            for x = 2:maxLevel
                %Build MS operators
                [CoarseGrid(x).MsR, CoarseGrid(x).MsP, CoarseGrid(x).C] = MSFVOperators(CoarseGrid(x-1), CoarseGrid(x), CoarseGrid(x-1).A_c, x);
                %Build coarse system (with MsFV)
                CoarseGrid(x).A_c = CoarseGrid(x).MsR*CoarseGrid(x-1).A_c*CoarseGrid(x).MsP;
                CoarseGrid(x).constant = 0;
            end
        end
        function BuildADMOperators(obj)
        end
    end
end