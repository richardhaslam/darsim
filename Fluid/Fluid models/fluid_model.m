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
            obj.Components = component();
        end
        function AddPhase(obj, Phase, index)
            obj.Phases(index) = Phase;
        end
        function AddComponent(obj, Comp, index)
            obj.Components(index) = Comp;
        end
        function Mob = ComputePhaseMobilities(obj, s)
            Mob = zeros(length(s), obj.NofPhases);
            kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, s);
            for i=1:obj.NofPhases
                Mob(:,i) = kr(:,i)/obj.Phases(i).mu;
            end
        end
        function z = ComputeMassFractions(s, x, rho)
            %Two phase, two component total mole fraction updater
            %Based on mass balance equation z_1 * rho_t = x11*rho1*s1 + x12*rho2*s2
            z(:,1) = (x(:,1).*rho(:,1).*s + x(:,2).*...
                rho(:,2).*(1-s))./(rho(:,1).*s + rho(:,2).*(1-s)); 
            z(:,2) = 1 - z(:,1);
        end
        function rhoT = ComputeTotalDensity(obj, s, rho)
            % Compute the total density
            rhoT = rho(:, 1) .* s + rho(:, 2) .* (1 - s);
        end
        function Status = ComputePhaseSaturation(obj, Status, SinglePhase)
            Status.S = Status.rho(:,2).*(Status.x1(:,2) - Status.z(:,1))./(Status.rho(:,1).*(Status.z(:,1)...
                - Status.x1(:,1)) + Status.rho(:,2).*(Status.x1(:,2) - Status.z(:,1)));
            Status.S(SinglePhase.onlyvapor == 1) = 1;
            Status.S(SinglePhase.onlyliquid == 1) = 0;
        end
    end
    methods (Abstract)
        obj = InitializeReservoir(obj);
        obj = InitializeInjectors(obj);
    end
end