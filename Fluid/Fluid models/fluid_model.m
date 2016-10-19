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
        function dMob = DMobDz(obj, Status, dSdz)
            dMobdS = obj.DMobDS(Status.S);
            % Use chain rule
            dMob(:,1) = dMobdS(:,1) .* dSdz;
            dMob(:,2) = dMobdS(:,2) .* dSdz;
        end
        function drho = DrhoDp(obj, p)
            drho = zeros(length(p), obj.NofPhases);
            for i=1:obj.NofPhases
                drho(:, i) = obj.Phases(i).DrhoDp(p);
            end
        end
        function drho = DrhoDz(obj, p, dxdz)
            drho = zeros(length(p), obj.NofPhases);
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
        function dPc = DPcDz(obj, Status)
            dPc = 0;
        end
        function dSdp = DSDp(obj, Status, drhodp, dni)
            ni = Status.ni;
            rhov = Status.rho(:,1);
            rhol = Status.rho(:,2);
            drhov = drhodp(:,1);
            drhol = drhodp(:,2);
            Num = rhol .* ni;
            Den = rhol .* ni + (1 - ni) .* rhov ;
            dNum = drhol .* ni + dni .* rhol;
            dDen = drhol .* ni + dni .* rhol + drhov - ni .* drhov - rhov .* dni;
            % final derivative
            dSdp = (dNum .* Den - dDen .* Num)./Den.^2;
        end
        function dSdz = DSDz(obj, Status, dni, dx1v, dx1l)
            rhov = Status.rho(:,1);
            rhol = Status.rho(:,2);
            x1v = Status.x1(:,1);
            x1l = Status.x1(:,2);
            ni = Status.ni;
            z = Status.z(:,1);
            % Derivative of S with respect to z
            Num = rhol .* (x1l - z);
            dNum = rhol .* (dx1l - 1);
            Den = rhov .* (z - x1v) + rhol .* (x1l - z); 
            dDen = rhov .* (1 - dx1v) + rhol .* (dx1l - 1);
            dSdz =(Den .* dNum - Num .* dDen) ./ Den.^2;
            Num1 = rhol .* ni;
            dNum1 = rhol .* dni;
            Den1 = rhol .* ni + (1 - ni) .* rhov;
            dDen1 = rhol .* dni - dni .* rhov;
            dSdz2 =(Den1 .* dNum1 - Num1 .* dDen1) ./ Den1.^2;
        end
        function drhotdz = DrhotDz(obj, Status, drho, dS)
            rho = Status.rho;
            S = Status.S;
            drhotdz = drho(:,1) .* S + rho(:,1) .* dS + drho(:,2) .* (1 - S) - rho(:,2) .* dS;
            % When it s one phase derivative is zero
            drhotdz(Status.S == 1) = 0;
            drhotdz(Status.S == 0) = 0;
        end
        function drhotdp = DrhotDp(obj, Status, drho, dS)
            rho = Status.rho;
            S = Status.S;
            drhotdp = drho(:,1) .* S + rho(:,1) .* dS + drho(:,2) .* (1 - S) - rho(:,2) .* dS;
            % When it s one phase derivative is zero
            drhotdp(Status.S == 1) = drho(Status.S == 1, 1);
            drhotdp(Status.S == 0) = drho(Status.S == 0, 2);
        end
    end
    methods (Abstract)
        obj = InitializeReservoir(obj);
        obj = InitializeInjectors(obj);
    end
end