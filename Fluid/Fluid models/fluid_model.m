% Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 17 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fluid_model < handle
    properties
        NofPhases
        NofComp
        Phases
        Components
        RelPermModel
    end
    methods
        function obj = fluid_model(n_phases, n_comp)
            obj.NofPhases = n_phases;
            obj.NofComp = n_comp;
            obj.Phases = phase();
        end
        function AddPhase(obj, Phase, index)
            obj.Phases(index) = Phase;
        end
        function Mob = ComputePhaseMobilities(obj, s)
            Mob = zeros(length(s), obj.NofPhases);
            kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, s);
            for i=1:obj.NofPhases
                Mob(:,i) = kr(:,i)/obj.Phases(i).mu;
            end
        end
    end
    methods (Abstract)
        obj = InitializeReservoir(obj);
        obj = InitializeInjectors(obj);
    end
end