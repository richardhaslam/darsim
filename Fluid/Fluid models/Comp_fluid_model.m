% Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 15 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Comp_fluid_model < fluid_model
    properties
        name = 'Compositional';
        FlashSettings
        KvaluesCalculator
    end
    methods
        function obj = Comp_fluid_model(n_phases, n_comp)
            obj@fluid_model(n_phases, n_comp);
        end
        function AddFlash(obj, flash)
            obj.FlashSettings = flash;
        end
        function InitializeReservoir(obj, Status)
            % Define initial values
            P_init = ones(length(Status.p), 1)*0.4e7;
            z_init = ones(length(Status.p), 1)*0.6;
            
            % 1. Assign initial valus
            Status.p = Status.p .* P_init;
            Status.z(:,1) = z_init;
            Status.z(:,2) = 1 - z_init;
            
            % 2. Update Composition of the phases (Flash)
            SinglePhase = obj.Flash(Status);
            
            % 2.a Compute Phase Density
            for i=1:obj.NofPhases
                Status.rho(:, i) = obj.Phases(i).ComputeDensity(Status.p);
            end
            
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
                SinglePhase = obj.Flash(Inj(i));
                for ph=1:obj.NofPhases
                    Inj(i).rho(:, ph)= obj.Phases(ph).ComputeDensity(Inj(i).p);
                end
                obj.ComputePhaseSaturation(Inj(i), SinglePhase);
                Inj(i).x2 = 1 - Inj(i).x1;
                Inj(i).Mob = obj.ComputePhaseMobilities(Inj(i).S);   
            end            
        end
        function SinglePhase = Flash(obj, Status)
            % Define SinglePhase objects
            SinglePhase.onlyliquid = zeros(length(Status.p), 1);
            SinglePhase.onlyvapor = zeros(length(Status.p), 1);
             
            k = obj.KvaluesCalculator.Compute(Status.p, Status.T, obj.Components);
            
            %% 2 Chek if we are in 2 phase region
            % 2.a: checking if it 's all liquid: checks if mix is below bubble
            % point
            % Transform Mass fractions to mol fractions
            z(:,1) =  Status.z(:,1) .* obj.Components(2).MM  ./ (obj.Components(2).MM * Status.z(:,1) +  obj.Components(1).MM * Status.z(:,2));
            z(:,2) =  Status.z(:,2) .* obj.Components(1).MM  ./ (obj.Components(2).MM * Status.z(:,1) +  obj.Components(1).MM * Status.z(:,2));
            %max(abs(z - Status.z))
            %z = Status.z;
            
            BubCheck = z .* k;
            BubCheck = sum(BubCheck, 2);
            
            x(BubCheck < 1, 2) = z (BubCheck < 1, 1);
            x(BubCheck < 1, 1) = 1;                     % This is to avoid having singular Jacobian matrix.
            SinglePhase.onlyliquid(BubCheck < 1) = 1;
            
            % 2.b: checking if it 's all vapor: checks if mix is above dew
            % point
            DewCheck = z ./ k;
            DewCheck = sum(DewCheck, 2);
            x(DewCheck < 1, 1) = z (DewCheck < 1, 1);
            x(DewCheck < 1, 2) = 1;                    % This is to avoid having singular Jacobian matrix.
            SinglePhase.onlyvapor(DewCheck < 1) = 1;
            
            %% 3. Actual Flash: solves for fv (vapor fraction)
            TwoPhase = ones(length(Status.p), 1);
            TwoPhase(SinglePhase.onlyliquid == 1) = 0;
            TwoPhase(SinglePhase.onlyvapor == 1) = 0;
            
            alpha = 0.5;
            
            %Initilaize variables
            fv = .5 * ones(length(Status.p),1);   % 50-50 split as inital guess
            %Single phase cells do not need to flash
            fv(SinglePhase.onlyvapor == 1) = 1;
            fv(SinglePhase.onlyliquid == 1)= 0;
            
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
                    if norm(h, inf) < obj.FlashSettings.TolFlash
                        
                        converged = 1;
                    end
                    itCounter = itCounter + 1;
                end
                alpha = alpha/2;
            end
            if ~converged
                disp('Warning: The flash calculator could not converge.');
            end
            
            %5. Solve for x's and y's
            % Solves for mole fractions in liquid phase
            x(TwoPhase == 1, 2) = z(TwoPhase == 1, 1) ./ (fv(TwoPhase == 1, 1) .* (k(TwoPhase == 1, 1) - 1) + 1);
            % Solves for mole fractions in gas phase
            x(TwoPhase == 1, 1) = k(TwoPhase == 1, 1) .* x(TwoPhase == 1, 2);
            
            % Convert x to mass fraction
            %x_old = x;
            x (TwoPhase == 1, 1) = x(TwoPhase == 1, 1) .* obj.Components(1).MM  ./ (obj.Components(1).MM * x(TwoPhase == 1,1) +  obj.Components(2).MM * (1 - x(TwoPhase == 1,1)));
            x (TwoPhase == 1, 2) = x(TwoPhase == 1, 2) .* obj.Components(1).MM  ./ (obj.Components(1).MM * x(TwoPhase == 1,2) +  obj.Components(2).MM * (1 - x(TwoPhase == 1,2)));
            %max(abs(x - x_old))
            
            % Copy it into Status object
            Status.x1 = x;
        end
        function SinglePhase = FindSinglePhaseCells(obj, Status, k)
            
            SinglePhase = zeros(length(Status.p), 1);
            
            BubCheck = Status.z .* k;
            BubCheck = sum(BubCheck, 2);
            
            Status.x1(BubCheck < 1, 2) = Status.z (BubCheck < 1, 1);
            Status.x1(BubCheck < 1, 1) = 1;                     % This is to avoid having singular Jacobian matrix.
            SinglePhase(BubCheck < 1) = 2;
            
            % 2.b: checking if it 's all vapor: checks if mix is above dew
            % point
            DewCheck = Status.z ./ k;
            DewCheck = sum(DewCheck, 2);
            Status.x1(DewCheck < 1, 1) = Status.z (DewCheck < 1, 1);
            Status.x1(DewCheck < 1, 2) = 1;                    % This is to avoid having singular Jacobian matrix.
            SinglePhase(DewCheck < 1) = 1;
        end
    end
end