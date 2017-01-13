% Compositional fluid model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 14 July 2016
%Last modified: 28 September 2016
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
            injvector = zeros(1, obj.NofComp);
            injvector(1) = 0.8;
            injvector(2) = 0.2;
            injvector(3) = 0.1;
            injxvector = [1 0 0 1 0 1];
            % Loop over all injectors
            for i=1:length(Inj)
                Inj(i).z = injvector;
                Inj(i).ni = 0.5;
                Inj(i).x = injxvector;
                SinglePhase = obj.Flash(Inj(i));                
                obj.ComputePhaseDensities(Inj(i));
                obj.ComputePhaseSaturation(Inj(i), SinglePhase);
                Inj(i).Mob = obj.ComputePhaseMobilities(Inj(i).S);   
            end            
        end
        function SinglePhase = Flash(obj, Status)
           SinglePhase = obj.FlashCalculator.Flash(Status, obj.Components, obj.Phases);
        end
        function ComputePhaseDensities(obj, Status)
            for i=1:obj.NofPhases
                Status.rho(:, i) = obj.Phases(i).ComputeDensity(Status, obj.Components);
            end
        end
        function SinglePhase = CheckNumberOfPhases(obj, SinglePhase, PreviousSinglePhase, Status, k)
            z = Status.z;
            % Check if single component
            SinglePhase(z(:, 1) == 1) = 1;
            SinglePhase(z(:, 1) == 0) = 2;
            Status.x(z(:, 1) == 1, 1) = 1;
            Status.x(z(:, 1) == 0, 1) = 1;
            Status.x(z(:, 1) == 1, 2) = 0;
            Status.x(z(:, 1) == 0, 2) = 0;
            Status.x(:,3:4) = 1 - Status.x(:,1:2);
            
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
        function dkdp = DKvalDp(obj, p)
            dkdp = obj.FlashCalculator.KvaluesCalculator.DKvalDp(p, obj.Components, obj.Phases);
        end
        function [dxdp, dxdz] = DxDpDz(obj, Status, SinglePhase)
            % x = x1v, x1l, x2v, x2l, ni
            % z = p, z1
            x = Status.x;
            k = obj.ComputeKvalues(Status);
            dk = obj.DKvalDp(Status);
            ni = Status.ni;
            N = length(Status.p);
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
%             dxdz(:,5,1) = -dxdz(:,5,1);
%             dxdz(:,6,1) = -dxdz(:,6,1);
%             
%             dxdz(:,1,2) = -dxdz(:,1,2);
%             dxdz(:,2,2) = -dxdz(:,2,2);
%             dxdz(:,3,2) = -dxdz(:,3,2);
%             dxdz(:,4,2) = -dxdz(:,4,2);
%             dxdz(:,5,2) = -dxdz(:,5,2);
%             dxdz(:,6,2) = -dxdz(:,6,2);
%             dxdz(:,7,2) = -dxdz(:,7,2);
            
            dxdz(SinglePhase == 1, 1, 1) = 1;
            dxdz(SinglePhase == 1, 2, 1) = 0;
            dxdz(SinglePhase == 1, 3, 1) = 0;
            dxdz(SinglePhase == 1, 4, 1) = 0;
            dxdz(SinglePhase == 1, 5, 1) = -1;
            dxdz(SinglePhase == 1, 6, 1) = 0;
            dxdz(SinglePhase == 1, 7, 1) = 0;
            
            dxdz(SinglePhase == 1, 1, 2) = 0;
            dxdz(SinglePhase == 1, 2, 2) = 0;
            dxdz(SinglePhase == 1, 3, 2) = 1;
            dxdz(SinglePhase == 1, 4, 2) = 0;
            dxdz(SinglePhase == 1, 5, 2) = -1;
            dxdz(SinglePhase == 1, 6, 1) = 0;
            dxdz(SinglePhase == 1, 7, 1) = 0;
            % 
            dxdz(SinglePhase == 2, 1, 1) = 0;
            dxdz(SinglePhase == 2, 2, 1) = 1;
            dxdz(SinglePhase == 2, 3, 1) = 0;
            dxdz(SinglePhase == 2, 4, 1) = 0;
            dxdz(SinglePhase == 2, 5, 1) = 0;
            dxdz(SinglePhase == 2, 6, 1) = -1;
            dxdz(SinglePhase == 2, 7, 1) = 0;
            
            dxdz(SinglePhase == 2, 1, 2) = 0;
            dxdz(SinglePhase == 2, 2, 2) = 0;
            dxdz(SinglePhase == 2, 3, 2) = 0;
            dxdz(SinglePhase == 2, 4, 2) = 1;
            dxdz(SinglePhase == 2, 5, 2) = 0;
            dxdz(SinglePhase == 2, 6, 2) = -1;
            dxdz(SinglePhase == 2, 7, 2) = 0;
            
        end
    end
end