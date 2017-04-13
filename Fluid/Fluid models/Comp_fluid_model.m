% Compositional fluid model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 14 July 2016
%Last modified: 7 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Comp_fluid_model < fluid_model
    properties
        FlashCalculator
    end
    methods
        function obj = Comp_fluid_model(n_phases, n_comp)
            obj@fluid_model(n_phases, n_comp);
            obj.name = 'Compositional';
        end
        function InitializeInjectors(obj, Inj)      
            % Loop over all injectors
            for i=1:length(Inj)
                state = status();
                % this is a temporary solution to a problem
                state.Properties('P_2') = property(length(Inj(i).Cells), 1, 0, 0, 0, Inj(i).p);
                state.Properties('S_1') = property(length(Inj(i).Cells), 1, 0,  0, 0, 1);
                state.Properties('S_2')= property(length(Inj(i).Cells), 1,  0, 0, 0, 1);
                state.Properties('z_1') = property(length(Inj(i).Cells), 1,  0, 0, 0, 1);
                state.Properties('z_2') = property(length(Inj(i).Cells), 1,  0, 0, 0, 0);
                state.Properties('rho_1') = property(length(Inj(i).Cells), 1,  0, 0, 0, 1);
                state.Properties('rho_2') = property(length(Inj(i).Cells), 1,  0, 0, 0, 0);
                state.Properties('x_1ph1') = property(length(Inj(i).Cells), 1,  0, 0, 0,1);
                state.Properties('x_1ph2') = property(length(Inj(i).Cells), 1,  0, 0, 0, 0);
                state.Properties('x_2ph1') = property(length(Inj(i).Cells), 1,  0, 0, 0, 0);
                state.Properties('x_2ph2') = property(length(Inj(i).Cells), 1,  0, 0, 0, 1);
                state.Properties('ni_1') = property(length(Inj(i).Cells), 1, 0,  0, 0, 0.5);
                state.Properties('ni_2') = property(length(Inj(i).Cells), 1, 0,  0, 0, 0.5);
                SinglePhase = obj.Flash(state);                
                obj.ComputePhaseDensities(state);
                obj.ComputePhaseSaturation(state, SinglePhase);
                Inj(i).Mob = obj.ComputePhaseMobilities(state.Properties('S_1').Value);
                % Copy values
                %Inj(i).z(:, 1) = state.Properties('z_1').Value;
                %Inj(i).z(:, 2) = state.Properties('z_2').Value;
                Inj(i).x(:, 1) = state.Properties('x_1ph1').Value;
                Inj(i).x(:, 2) = state.Properties('x_2ph1').Value;
                Inj(i).x(:, 3:4) = 1 - Inj(i).x(:,1:2);
                Inj(i).rho(:, 1) = state.Properties('rho_1').Value;
                Inj(i).rho(:, 2) = state.Properties('rho_2').Value;
            end            
        end
        function SinglePhase = Flash(obj, Status)
           SinglePhase = obj.FlashCalculator.Flash(Status, obj.Components, obj.Phases);
        end
        function ComputePhaseSaturation(obj, Status, SinglePhase)
            N = length(SinglePhase);
            rho = zeros(N, obj.NofPhases);
            ni = zeros(N, obj.NofPhases);
            % Store in local variables
            NiRhoT = zeros(N, 1);
            for i=1:obj.NofPhases
                rho(:,i) = Status.Properties(['rho_',num2str(i)]).Value;
                ni(:,i) = Status.Properties(['ni_',num2str(i)]).Value;
                NiRhoT = NiRhoT + ni(:,i)./rho(:,i);
            end
            % Compute saturations
            for i=1:obj.NofPhases
                S = Status.Properties(['S_',num2str(i)]);
                S.Value = (ni(:,i)./rho(:,i)) ./ NiRhoT;
            end
        end
        function ComputePhaseDensities(obj, Status, SinglePhase)
            for i=1:obj.NofPhases
                rho = Status.Properties(strcat('rho_', num2str(i)));
                rho.Value = obj.Phases(i).ComputeDensity(Status.Properties('P_2').Value, obj.Components);
            end
        end
        function z = ComputeTotalFractions(obj, Status, N)
            %Two phase, two component total mole fraction updater
            %Based on mass balance equation z_1 * rho_t = x11*rho1*s1 + x12*rho2*s2
            S = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            x = zeros(N, obj.NofComp*obj.NofPhases);
            
            for ph=1:obj.NofPhases
                rho(:, ph) = Status.Properties(['rho_', num2str(ph)]).Value;
                S(:, ph) = Status.Properties(['S_', num2str(ph)]).Value;
                for c=1:obj.NofComp
                    x(:,(c-1)*obj.NofPhases + ph) = Status.Properties(['x_', num2str(c),'ph',num2str(ph)]).Value;
                end
            end
            rhoT = Status.Properties('rhoT').Value;
            for i=1:obj.NofComp
                z = Status.Properties(['z_', num2str(i)]);
                Num = zeros(N,1);
                for j=1:obj.NofPhases
                    Num = Num + x(:,(i-1)*obj.NofPhases + j) .* S(:, j) .* rho(:,j); 
                end
                z.Value = Num ./ rhoT; 
            end
        end
        function SinglePhase = CheckNumberOfPhases(obj, SinglePhase, PreviousSinglePhase, Status, k)
            N = length(SinglePhase);
            z = zeros(N, obj.NofComp);
            for i=1:obj.NofComp
                z(:, i) = Status.Properties(['z_', num2str(i)]).Value;
                for j=1:obj.NofPhases
                    x = Status.Properties(['x_', num2str(i),'ph',num2str(j)]);
                    x.Value(z(:,i) == 1) = i==j;
                    x.Value(z(:,i) == 0) = i==j;
                end
            end
            
            % Check if single component
            SinglePhase(z(:,1) == 1) = 1;
            SinglePhase(z(:,1) == 0) = 2;
            
            % Transform mass fractions into mole fractions to check phase
            % state
            BubCheck = zeros(length(z), 2);
            BubCheck(PreviousSinglePhase == 2, :) = z(PreviousSinglePhase == 2,:) .* k (PreviousSinglePhase == 2, :);
            BubCheck = sum(BubCheck, 2);
            SinglePhase(BubCheck > 1) = 0;
            
            % 2.b: checking if it 's all vapor: checks if mix is above dew
            % point
            DewCheck = zeros(length(z), 2);
            DewCheck(PreviousSinglePhase == 1,:) = z(PreviousSinglePhase == 1,:) ./ k(PreviousSinglePhase == 1,:);
            DewCheck = sum(DewCheck, 2);
            SinglePhase(DewCheck(PreviousSinglePhase == 1) > 1) = 0;
        end
        function k = ComputeKvalues(obj, Status)
            k = obj.FlashCalculator.KvaluesCalculator.Compute(Status, obj.Components, obj.Phases);
        end
        function dkdp = DKvalDp(obj, Status)
            dkdp = obj.FlashCalculator.KvaluesCalculator.DKvalDp(Status, obj.Components, obj.Phases);
        end
        function [dxdp, dxdz] = DxDpDz(obj, Status, SinglePhase)
            % x = x1v, x1l, x2v, x2l, ni
            % z = p, z1
            N = length(SinglePhase);
            Ni = zeros(N, obj.NofPhases);
            x = zeros(N, obj.NofPhases*obj.NofComp);
            for j=1:obj.NofPhases
                Ni(:, j) = Status.Properties(['ni_', num2str(j)]).Value;
                for i=1:obj.NofComp
                    x(:,(i-1)*obj.NofPhases + j) = Status.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
                end
            end
            ni = Ni(:,1);
            
            k = obj.ComputeKvalues(Status);
            dk = obj.DKvalDp(Status);
            % Loop over all cells and do local inversion
            dxdp = zeros(N, 2*obj.NofComp+1);
            dxdz = zeros(N, 2*obj.NofComp+1, obj.NofComp-1);
            for i = 1:N
                dFdx = zeros(2*obj.NofComp + 1);
                dFdz = zeros(2*obj.NofComp + 1, obj.NofComp);
                for c=1:obj.NofComp
                    % x_cv - k_c x_cl = 0
                    dFdx (c, (c-1)*2+1:(c-1)*2+2) = [1, -k(i,c)];
                    dFdz (c, 1) = -dk(i,c) * x(i, (c-1)*2+2);
                    % z_c - ni*x_cv - (1-ni)*x_cl = 0
                    dFdx (obj.NofComp + c, (c-1)*2+1:(c-1)*2+2) = [-ni(i), ni(i)-1];
                    dFdx (obj.NofComp + c, end) = x(i, (c-1)*2+2) - x(i, (c-1)*2+1);
                    if c<obj.NofComp
                        dFdz (obj.NofComp + c, c+1) = 1;
                    end
                    % Sum x_cv-x_cl = 0
                    dFdx (2*obj.NofComp + 1, (c-1)*2+1:(c-1)*2+2) = [1, -1]; 
                end
                dFdz (2*obj.NofComp, 2:obj.NofComp) = -1;
                dxdp(i,:) = (dFdx\dFdz(:,1))';
                for j=1:obj.NofComp-1
                    dxdz(i,:,j) = (dFdx\dFdz(:,j+1))';
                end
            end
            dxdp(SinglePhase == 1, :) = 0;
            dxdp(SinglePhase == 2, :) = 0;
            
            dxdz(SinglePhase == 1, 1) = 1; %1
            dxdz(SinglePhase == 1, 2) = 0;
            dxdz(SinglePhase == 1, 3) = -1; %-1
            dxdz(SinglePhase == 1, 4) = 0;
            dxdz(SinglePhase == 1, 5) = 0;
            
            % 
            dxdz(SinglePhase == 2, 1) = 0;
            dxdz(SinglePhase == 2, 2) = 1;%1
            dxdz(SinglePhase == 2, 3) = 0;
            dxdz(SinglePhase == 2, 4) = -1;%-1
            dxdz(SinglePhase == 2, 5) = 0;
            
        end
        function drho = DrhoDz(obj, Status, SinglePhase)
            N = length(Status.Properties('P_2').Value);
            drho = zeros(N, obj.NofPhases);
        end
        function dSdp = DSDp(obj, Status, drhodp, dni)
            ni = Status.Properties('ni_1').Value;
            rhov = Status.Properties('rho_1').Value;
            rhol = Status.Properties('rho_2').Value;
            drhov = drhodp(:,1);
            drhol = drhodp(:,2);
            Num = rhol .* ni;
            Den = rhol .* ni + (1 - ni) .* rhov ;
            dNum = drhol .* ni + dni .* rhol;
            dDen = drhol .* ni + dni .* rhol + drhov - ni .* drhov - rhov .* dni;
            % final derivative
            dSdp = (dNum .* Den - dDen .* Num)./Den.^2;
        end
        function dMobdp = DMobDp(obj, Status, dSdp)
            dMobdS = obj.DMobDS(Status.Properties('S_1').Value); 
            dMobdp = dMobdS .* dSdp;
        end
        function dSdz = DSDz(obj, Status, dni)
            % Derivative of S with respect to overall composition
            rhov = Status.Properties('rho_1').Value;
            rhol = Status.Properties('rho_2').Value;
            ni = Status.Properties('ni_1').Value;
            dSdz = zeros(length(ni), obj.NofComp - 1);
            for j = 1:obj.NofComp-1
                Num1 = rhol .* ni;
                dNum1 = rhol .* dni(:,j);
                Den1 = rhol .* ni + (1 - ni) .* rhov;
                dDen1 = rhol .* dni(:,j) - dni(:,j) .* rhov;
                dSdz(:, j) =(Den1 .* dNum1 - Num1 .* dDen1) ./ Den1.^2;
            end
            dSdz(Status.Properties('S_1').Value == 0, 2) = 0;
            dSdz(Status.Properties('S_1').Value == 0, 1) = 0;
            dSdz(Status.Properties('S_1').Value == 1, 2) = 0;
            dSdz(Status.Properties('S_1').Value == 1, 1) = 0;
        end
        function drhotdz = DrhotDz(obj, Status, drho, dS)
            N = length(Status.Properties('rho_1').Value);
            rho = zeros(N, obj.NofPhases);
            S = zeros(N, obj.NofPhases);
            for i=1:obj.NofPhases
                rho(:,i) = Status.Properties(['rho_',num2str(i)]).Value;
                S(:,i) = Status.Properties(['S_',num2str(i)]).Value;
            end
            drhotdz = zeros(length(S), obj.NofComp-1);
            for j=1:obj.NofComp-1
                drhotdz(:,j) = drho(:,1) .* S(:,1) + rho(:,1) .* dS(:,j) + drho(:,2) .* S(:,2) - rho(:,2) .* dS(:,j);
            end
            % When it s one phase
            drhotdz(S(:,1) == 1, :) = 0;
        end
        function dPc = DPcDz(obj, Status, dSdz)
            dPcdS = obj.DPcDS(Status.Properties('S_1').Value);
            dPc = zeros(length(dPcdS), obj.NofComp-1);
            for j=1:obj.NofComp-1
                % Use chain rule
                dPc(:,j) = dPcdS(:,1) .* dSdz(:,j);
            end
        end
    end
end