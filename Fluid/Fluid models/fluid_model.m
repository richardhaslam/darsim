% Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 19 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fluid_model < handle
    properties
        NofPhases
        NofComp
        Phases
        Components
        RelPermModel
        CapillaryModel
        WettingPhaseIndex
        GravityModel
    end
    methods
        function obj = fluid_model(n_phases, n_comp)
            obj.NofPhases = n_phases;
            obj.NofComp = n_comp;
            obj.Phases = phase();
            obj.Components = component();
            obj.WettingPhaseIndex = 2;
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
        function z = ComputeTotalFractions(obj, s, x, rho)
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
        function ComputePhaseSaturation(obj, Status, SinglePhase)
            Status.S = Status.rho(:,2).*(Status.x1(:,2) - Status.z(:,1))./(Status.rho(:,1).*(Status.z(:,1)...
                - Status.x1(:,1)) + Status.rho(:,2).*(Status.x1(:,2) - Status.z(:,1)));
            Status.S(SinglePhase.onlyvapor == 1) = 1;
            Status.S(SinglePhase.onlyliquid == 1) = 0;
        end
        function dMob = DMobDS(obj, s)
            dMob = zeros(length(s), obj.NofPhases);
            dkr = obj.RelPermModel.ComputeDerivative(obj.Phases, s);
            for i=1:obj.NofPhases
                dMob(:,i) = dkr(:,i)/obj.Phases(i).mu;
            end
        end
        function drho = DrhoDp(obj, p)
            drho = zeros(length(p), obj.NofPhases);
            for i=1:obj.NofPhases
                drho(:, i) = obj.Phases(i).DrhoDp(p);
            end
        end
        function Pc = ComputePc(obj, S)
            switch(obj.WettingPhaseIndex)
                case(1)
                    S = (S - obj.Phases(2).sr)./(1 - obj.Phases(2).sr) + 0.1;
                    Pc = obj.CapillaryModel.ComputePc(S);
                case(2)
                    S = 1 - S;
                    S = (S - obj.Phases(1).sr)./(1 - obj.Phases(1).sr) + 0.1;
                    Pc = -obj.CapillaryModel.ComputePc(S);
            end   
        end
        function dPc = DPcDS(obj, S)
            switch(obj.WettingPhaseIndex)
                case(1)
                    S = (S - obj.Phases(2).sr)./(1 - obj.Phases(2).sr) + 0.1;
                    dPc = obj.CapillaryModel.dPcdS(S);
                case(2)
                    S = 1 - S;
                    S = (S - obj.Phases(1).sr)./(1 - obj.Phases(1).sr) + 0.1;
                    dPc =  - obj.CapillaryModel.dPcdS(S);
            end
            
        end
    end
    methods (Abstract)
        obj = InitializeReservoir(obj);
        obj = InitializeInjectors(obj);
    end
end