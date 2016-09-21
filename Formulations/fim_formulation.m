%  FIM Formulation base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 19 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fim_formulation < handle
    properties
        NumOfEquations
        UpWindPh1
        UpWindPh2
        Uph1
        Uph2
        Tph1
        Tph2
        Mob
        dMob
        dPc
        drho
    end
    methods (Abstract)
        obj = BuildResidual(obj)
        obj = BuildJacobian(obj)
        obj = UpdateState(obj)
        obj = Reset(obj)
    end
    methods
        function UpWindAndPhaseRockFluxes(obj, Grid, Phases, Pot)
            [obj.UpWindPh1, obj.Uph1] = Phases(1).UpWindAndRockFluxes(Grid, Pot(:,1));
            [obj.UpWindPh2, obj.Uph2] = Phases(2).UpWindAndRockFluxes(Grid, Pot(:,2));
        end
    end
end