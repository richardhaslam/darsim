% Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 16 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fluid_model < handle
    properties
        name
        NofPhases
        NofComp
        Phases
        Components
        RelPermModel
        CapillaryModel
        WettingPhaseIndex
    end
    methods
        function obj = fluid_model(n_phases, n_comp)
            obj.NofPhases = n_phases;
            obj.NofComp = n_comp;
            obj.Phases = phase.empty;
            obj.Components = component.empty;
            obj.WettingPhaseIndex = 1;
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
        function ComputeTotalDensity(obj, State)
            % Compute the total density
            s = State.Properties('S_1');
            rho = State.Properties('rho_1');
            rhoT = rho.Value .* s.Value;
            for i=2:obj.NofPhases
                s = State.Properties(['S_', num2str(i)]);
                rho = State.Properties(['rho_', num2str(i)]);
                rhoT = rhoT + rho.Value .* s.Value;
            end
            RhoT = State.Properties('rhoT');
            RhoT.Value = rhoT;
        end
        function dMob = DMobDS(obj, s)
            dMob = zeros(length(s), obj.NofPhases);
            dkr = obj.RelPermModel.ComputeDerivative(obj.Phases, s);
            for i=1:obj.NofPhases
                dMob(:,i) = dkr(:,i)/obj.Phases(i).mu;
            end
        end
        function dMob = DMobDz(obj, Status, dSdz)
            dMobdS = obj.DMobDS(Status.Properties('S_1').Value);
            dMob = zeros(length(dMobdS), obj.NofComp-1);
            for j=1:obj.NofComp-1
                % Use chain rule
                dMob(:,1,j) = dMobdS(:,1) .* dSdz(:,j);
                dMob(:,2,j) = dMobdS(:,2) .* dSdz(:,j);
            end
        end
        function drho = DrhoDp(obj, Status, SinglePhase)
            drho = zeros(length(Status.Properties('P_1').Value), obj.NofPhases);
            for i=1:obj.NofPhases
                drho(:, i) = obj.Phases(i).DrhoDp(Status.Properties(strcat('P_', num2str(obj.NofPhases))).Value);
            end
        end
        function ComputePc(obj, Status)
            % Define local handles
            S = Status.Properties('S_1').Value;
            P1 = Status.Properties('P_1');
            P2 = Status.Properties('P_2');
            Pc = Status.Properties('Pc');
            % Compute Pc as a function of S_1
            switch(obj.WettingPhaseIndex)
                case(1)
                    S = (S - obj.Phases(1).sr)./(1 - obj.Phases(1).sr) + 0.1;
                    S (S < obj.Phases(1).sr) = 0.1;
                    Pc.Value = obj.CapillaryModel.ComputePc(S);
                case(2)
                    S = 1 - S;
                    S = (S - obj.Phases(2).sr)./(1 - obj.Phases(2).sr) + 0.1;
                    Pc.Value = -obj.CapillaryModel.ComputePc(S);
            end
            % Update P_1
            P1.Value = P2.Value - Pc.Value;
        end
        function dPc = DPcDS(obj, S)
            switch(obj.WettingPhaseIndex)
                case(1)
                    S = (S - obj.Phases(1).sr)./(1 - obj.Phases(1).sr) + 0.1;
                    dPc = obj.CapillaryModel.dPcdS(S);
                    dPc (S < obj.Phases(1).sr) = 0.0;
                case(2)
                    S = 1 - S;
                    S = (S - obj.Phases(2).sr)./(1 - obj.Phases(2).sr) + 0.1;
                    dPc =  - obj.CapillaryModel.dPcdS(S);
            end
            
        end
        function drhotdp = DrhotDp(obj, Status, drho, dS)
            N = length(Status.Properties('rho_1').Value);
            rho = zeros(N, obj.NofPhases);
            S = zeros(N, obj.NofPhases);
            for i=1:obj.NofPhases
                rho(:,i) = Status.Properties(['rho_',num2str(i)]).Value;
                S(:,i) = Status.Properties(['S_',num2str(i)]).Value;
            end
            drhotdp = drho(:,1) .* S(:,1) + rho(:,1) .* dS + drho(:,2) .* S(:,2) - rho(:,2) .* dS;
            % When it s one phase derivative is the one of the existing
            % phase
            drhotdp(S(:,1) == 1) = drho(S(:,1) == 1, 1);
            drhotdp(S(:,1) == 0) = drho(S(:,1) == 0, 2);
        end
    end
    methods (Abstract)
        obj = Flash(obj);
        obj = InitializeInjectors(obj);
    end
end