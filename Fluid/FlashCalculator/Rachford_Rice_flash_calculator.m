% Rachford-Rice Flash Calculator base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 26 October 2016
%Last modified: 26 October 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Rachford_Rice_flash_calculator < Kvalues_flash_calculator
    properties
        tol = 1e-5;
    end
    methods
        function SinglePhase = Flash(obj, Status, Components, Phases)
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
            
            %% 2 Compute K-values
            k = obj.KvaluesCalculator.Compute(Status, Components, Phases);
            
            %% 3 Chek if we are in 2 phase region
            % 2.a: checking if it 's all liquid: checks if mix is below bubble
            BubCheck = z .* k;
            BubCheck = sum(BubCheck, 2);
            
            x(BubCheck - 1 <= obj.tol, 2) = z (BubCheck - 1 < obj.tol, 1);
            x(BubCheck - 1 <= obj.tol, 1) = 1;                     % This is to avoid having singular Jacobian matrix.
            SinglePhase(BubCheck - 1 <= obj.tol) = 2; % It s all liquid
            
            % 2.b: checking if it 's all vapor: checks if mix is above dew
            % point
            DewCheck = z ./ k;
            DewCheck = sum(DewCheck, 2);
            x(DewCheck - 1 < obj.tol, 1) = z (DewCheck - 1 < obj.tol, 1);
            x(DewCheck - 1 < obj.tol, 2) = 1;                    % This is to avoid having singular Jacobian matrix.
            SinglePhase(DewCheck - 1 < obj.tol) = 1;  % It s all vapour
            
            %% 3. Actual Flash: solves for fv (vapor fraction)
            TwoPhase = ones(length(Status.p), 1);
            TwoPhase(SinglePhase > 0 ) = 0;
            
            
            %Initilaize variables
            alpha = ones(length(Status.p),1);
            fv = Status.ni;
            
            %Single phase cells do not need to flash
            fv(SinglePhase == 1) = 1;
            fv(SinglePhase == 2)= 0;
            
            % Find fv with the tangent method
            converged = 0;
            itLimit = 10000;
            while ~converged && min(alpha) > 0.01
                itCounter = 0;
                while itCounter < itLimit && ~converged
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
                    dh(TwoPhase == 0) = 1;
                    fvnew = alpha .* (-h ./ dh) + fv;
                    
                    fv = fvnew;
                    if norm(h, inf) < 1e-10
                        converged = 1;
                        disp(['Rachford-Rice converged in ', num2str(itCounter + 1), ' iterations, with alpha ', num2str(min(alpha))])
                    end
                    itCounter = itCounter + 1;
                end
                alpha (abs(h) > 1e-10) = alpha (abs(h) > 1e-10)/2;
                fv (fv > 1) = 0.9;
                fv (fv < 0) = 0.1;
            end
            if ~converged
                [~, cellIndex] = max(abs(h));
                disp('Warning: Flash did not fully converge!');
                disp(['The residual norm of the equilibrium equation is ', num2str(norm(h, inf)), ' in cell ', num2str(cellIndex)]);
                fv (fv > 1) = 1;
                fv(fv < 0 ) = 0;
            end
            
            %5. Solve for x's and y's
            % Solves for mole fractions in liquid phase
            x(TwoPhase == 1, 2) = z(TwoPhase == 1, 1) ./ (fv(TwoPhase == 1, 1) .* (k(TwoPhase == 1, 1) - 1) + 1);
            % Solves for mole fractions in gas phase
            x(TwoPhase == 1, 1) = k(TwoPhase == 1, 1) .* x(TwoPhase == 1, 2);
            
            % Copy it into Status object
            Status.x(:,1:2) = x;
            Status.x(:,3:4) = 1 - x;
            Status.ni = fv;
        end
    end
end