%  FIM Formulation base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 19 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fim_formulation < handle
    properties
        NumOfEquations
        UpWindW
        UpWindNw
        Uw
        Unw
    end
    methods (Abstract)
        obj = BuildResidual(obj)
        obj = BuildJacobian(obj)
        obj = UpdateState(obj)
    end
    methods
        function UpWindAndPhaseRockFluxes(obj, Grid, Phases, p, pc)
            [obj.UpWindW, obj.Uw] = Phases(1).UpWindAndRockFluxes(Grid, p - pc);
            [obj.UpWindNw, obj.Unw] = Phases(2).UpWindAndRockFluxes(Grid, p);
        end
    end
end