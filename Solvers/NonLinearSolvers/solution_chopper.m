% NL solver base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef solution_chopper < handle
    properties
        MaxRatio = 1;
        MaxDelta
    end
    methods
        function DefineMaxDelta(obj, x)
            obj.MaxDelta = obj.MaxRatio * x;
            obj.MaxDelta(obj.MaxDelta == 0) = .1;
        end
        function delta = Chop(obj, delta)
            %ratio = max(delta ./ obj.MaxDelta);
%             if ratio>1
%                 delta = 0.2*delta;
%             end
        end
    end
end