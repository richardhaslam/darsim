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
        KvaluesCalculator
    end
    methods
        function obj = Comp_fluid_model(n_phases, n_comp)
            obj@fluid_model(n_phases, n_comp);
            obj.name = 'Compositional';
        end
        function SinglePhase = InitializeReservoir(obj, Status)
            % Define initial values
            P_init = ones(length(Status.p), 1)*1e5;
            z_init = ones(length(Status.p), 1)*0.03738;
            
            % 1. Assign initial valus
            Status.p = Status.p .* P_init;
            Status.z(:,1) = z_init;
            Status.z(:,2) = 1 - z_init;
            
            % 2. Update Composition of the phases (Flash)
            k = obj.ComputeKvalues(Status);
            SinglePhase = obj.Flash(Status, k);
            
            % 2.a Compute Phase Density
            obj.ComputePhaseDensities(Status);
            
            %3. Update S based on components mass balance
            obj.ComputePhaseSaturation(Status, SinglePhase);
            
            % 4. Total Density
            Status.rhoT = obj.ComputeTotalDensity(Status.S, Status.rho);
            
            % 5. Compute initial Pc
            Status.pc = obj.ComputePc(Status.S);
        end
        function InitializeInjectors(obj, Inj)
            % Loop over all injectors
            for i=1:length(Inj)
                Inj(i).z = [1 0];
                k = obj.ComputeKvalues(Inj(i));
                SinglePhase = obj.Flash(Inj(i), k);                
                obj.ComputePhaseDensities(Inj(i));
                obj.ComputePhaseSaturation(Inj(i), SinglePhase);
                Inj(i).x2 = 1 - Inj(i).x1;
                Inj(i).Mob = obj.ComputePhaseMobilities(Inj(i).S);   
            end            
        end
        function SinglePhase = Flash(obj, Status, k)
            % Define SinglePhase objects
            SinglePhase = zeros(length(Status.p), 1);
            z = Status.z;
            
            %% 1 Check if there are 2 components
            x(z(:, 1) == 1, 1) = 1;  
            x(z(:, 1) == 1, 2) = 1;
            SinglePhase(z(:, 1) == 1) = 1; % All vapour
            x(z(:, 1) == 0, 1) = 1;  
            x(z(:, 1) == 0, 2) = 1;
            SinglePhase(z(:, 1) == 0) = 2; % All liquid
            
            %% 2 Chek if we are in 2 phase region
            % 2.a: checking if it 's all liquid: checks if mix is below bubble
            BubCheck = z .* k;
            BubCheck = sum(BubCheck, 2);
            
            x(BubCheck < 1, 2) = z (BubCheck < 1, 1);
            x(BubCheck < 1, 1) = 1;                     % This is to avoid having singular Jacobian matrix.
            SinglePhase(BubCheck < 1) = 2; % It s all liquid
            
            % 2.b: checking if it 's all vapor: checks if mix is above dew
            % point
            DewCheck = z ./ k;
            DewCheck = sum(DewCheck, 2);
            x(DewCheck < 1, 1) = z (DewCheck < 1, 1);
            x(DewCheck < 1, 2) = 1;                    % This is to avoid having singular Jacobian matrix.
            SinglePhase(DewCheck < 1) = 1;  % It s all vapour
            
            %% 3. Actual Flash: solves for fv (vapor fraction)
            TwoPhase = ones(length(Status.p), 1);
            TwoPhase(SinglePhase > 0 ) = 0;
            
            alpha = 0.5;
            
            %Initilaize variables
            fv = .5 * ones(length(Status.p),1);   % 50-50 split as inital guess
            %Single phase cells do not need to flash
            fv(SinglePhase == 1) = 1;
            fv(SinglePhase == 2)= 0;
            
            % Find fv with the tangent method
            converged = 0;
            while ~converged && alpha > 0.1
                itCounter = 0;
                while itCounter < 200 && ~converged
                    %Finds hi for each component
                    hi(:,1) = (z(:,1) .* k(:,1)) ./ (fv .* (k(:,1) - 1) + 1);
                    hi(:,2) = (z(:,2) .* k(:,2)) ./ (fv .* (k(:,2) - 1) + 1);
                    %Finds the derivative of hi for each component
                    dhi(:, 1) = (z(:,1) .* (k(:,1) - 1).^2) ./ ((fv .* (k(:,1) - 1) + 1).^2);
                    dhi(:, 2) = (z(:,2) .* (k(:,2) - 1).^2) ./ ((fv .* (k(:,2) - 1) + 1).^2);
                    h = sum(hi, 2) - 1;
                    dh = - sum(dhi, 2);
                    
                    %Update fv
                    h(TwoPhase == 0) = 0;
                    fvnew = alpha * (-h ./ dh) + fv;
                    
                    fv = fvnew;
                    if norm(h, inf) < 1e-10
                        converged = 1;
                    end
                    itCounter = itCounter + 1;
                end
                alpha = alpha/2;
            end
            if ~converged
                disp('Warning: Flash did not converge!');
                disp(['The residual norm of the equilibrium equation is ', num2str(norm(h, inf))]);
            end
            
            %5. Solve for x's and y's
            % Solves for mole fractions in liquid phase
            x(TwoPhase == 1, 2) = z(TwoPhase == 1, 1) ./ (fv(TwoPhase == 1, 1) .* (k(TwoPhase == 1, 1) - 1) + 1);
            % Solves for mole fractions in gas phase
            x(TwoPhase == 1, 1) = k(TwoPhase == 1, 1) .* x(TwoPhase == 1, 2);
            
            % Convert x to mass fraction
            %x_old = x;
            %x (TwoPhase == 1, 1) = x(TwoPhase == 1, 1)  ./ (x(TwoPhase == 1,1) +  (1 - x(TwoPhase == 1,1)));
            %x (TwoPhase == 1, 2) = x(TwoPhase == 1, 2) ./ (x(TwoPhase == 1,2) +  (1 - x(TwoPhase == 1,2)));
            %max(abs(x - x_old))
            
            % Copy it into Status object
            %% 1 Check if there are 2 components
            Status.x1 = x;
            Status.ni = fv .* sum(z, 2);
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
            Status.x1(z(:, 1) == 1, 1) = 1;
            Status.x1(z(:, 1) == 0, 1) = 1;
            Status.x1(z(:, 1) == 1, 2) = 0;
            Status.x1(z(:, 1) == 0, 2) = 0;
            
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
            k = obj.KvaluesCalculator.Compute(Status.p, Status.T, obj.Components);
        end
        function dkdp = DKvalDp(obj, p)
            dkdp = obj.KvaluesCalculator.DKvalDp(p);
        end
        function [dxdp, dxdz] = DxDpDz(obj, Status, SinglePhase)
            % Find sensitivities
            % x = x1v, x1l, x2v, x2l, ni
            % z = p, z1
            k = obj.ComputeKvalues(Status);
            dk = obj.DKvalDp(Status.p);
            ni = Status.ni;
            x1 = Status.x1(:,1);
            y1 = Status.x1(:,2);
            x2 = 1 - x1;
            y2 = 1 - y1;
            k1 = k(:,1);
            k2 = k(:,2);
            dk1 = dk(:,1);
            dk2 = dk(:,2);
            
            N = length(Status.p);
            % Loop over all cells and do inversion local inversion
            dxdp = zeros(N, 5);
            dxdz = zeros(N, 5);
            for i = 1:N
                dFdz = zeros(5,2);
                dFdx = zeros(5,5);
                % Equation 1
                % dF1/dx1v
                dFdx(1, 1) = 1;
                % dF1/dx1l
                dFdx(1, 2) = -k1(i);
                % Equation 2
                % dF2/dx2v
                dFdx(2, 3) = 1;
                % dF2/dx2l
                dFdx(2, 4) = -k2(i);
                % Equation 3
                % dF3/dx1v
                dFdx(3, 1) = - ni(i);
                % dF3/dx1l
                dFdx(3, 2) = -(1 - ni(i));
                % dF3dni
                dFdx(3, 5) = y1(i) - x1(i);
                % Equation 4
                % dF4/dx2v
                dFdx(4, 3) = - ni(i);
                % dF4/dx2l
                dFdx(4, 4) = - (1 - ni(i));
                % dF4dni
                dFdx(4, 5) = y2(i) - x2(i);
                % Equation 5
                % dF5/dx1v
                dFdx(5, 1) = 1;
                % dF5/dx1l
                dFdx(5, 2) = -1;
                % dF5/dx2v
                dFdx(5, 3) = 1;
                % dF5/dx2l
                dFdx(5, 4) = -1;
                % dFdz
                dFdz(1, 1) = -dk1(i) * y1(i);
                dFdz(2, 1) = -dk2(i) * y2(i);
                dFdz(3, 2) =  1;
                dFdz(4, 2) = -1;
                dFdx = sparse(dFdx);
                dxdp(i,:) = (dFdx\dFdz(:,1))';
                dxdz(i,:) = (dFdx\dFdz(:,2))';
            end
            dxdz(SinglePhase == 1, 1) = 1;
            dxdz(SinglePhase == 1, 2) = 0;
            dxdz(SinglePhase == 1, 3) = -1;
            dxdz(SinglePhase == 1, 4) = 0;
            dxdz(SinglePhase == 1, 5) = 0;
            
            % 
            dxdz(SinglePhase == 2, 1) = 0;
            dxdz(SinglePhase == 2, 2) = 1;
            dxdz(SinglePhase == 2, 3) = 0;
            dxdz(SinglePhase == 2, 4) = -1;
            dxdz(SinglePhase == 2, 5) = 0;
            
        end
    end
end