% Rachford-Rice Flash Calculator base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 26 October 2016
%Last modified: 6 APril 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Rachford_Rice_flash_calculator < Kvalues_flash_calculator
    properties
        tol = 1e-5;
    end
    methods
        function SinglePhase = Flash(obj, Status, Components, Phases)
            % Define SinglePhase objects
            nc = length(Components);
            nph = length(Phases);
            N = length(Status.Properties('z_1').Value);
            z = zeros(N, nc);
            SinglePhase = zeros(N, 1);
            for i=1:nc
                z(:, i) = Status.Properties(['z_', num2str(i)]).Value;
            end
            
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
            % 3.a: checking if it 's all liquid: checks if mix is below bubble
            BubCheck = z .* k;
            BubCheck = sum(BubCheck, 2);
            
            x(BubCheck - 1 <= obj.tol, 2) = z (BubCheck - 1 < obj.tol, 1);
            x(BubCheck - 1 <= obj.tol, 1) = 1; % This is to avoid having singular Jacobian matrix.
            SinglePhase(BubCheck - 1 <= obj.tol) = 2; % It s all liquid
            
            % 3.b: checking if it 's all vapor: checks if mix is above dew
            % point
            DewCheck = z ./ k;
            DewCheck = sum(DewCheck, 2);
            x(DewCheck - 1 < obj.tol, 1) = z (DewCheck - 1 < obj.tol, 1);
            x(DewCheck - 1 < obj.tol, 2) = 1; % This is to avoid having singular Jacobian matrix.
            SinglePhase(DewCheck - 1 < obj.tol) = 1;  % It s all vapour
            
            %% 4. Actual Flash: solves for fv (vapor fraction)
            TwoPhase = ones(length(SinglePhase), 1);
            TwoPhase(SinglePhase > 0 ) = 0;
            
            
            % Initilaize variables
            alpha = ones(N, 1);
            fv = Status.Properties('ni_1').Value;
            
            % Single phase cells do not need to flash
            fv(SinglePhase == 1) = 1;
            fv(SinglePhase == 2) = 0;
            
            % Find fv with the tangent method
%             converged = 0;
%             itLimit = 10000;
%             hi = zeros(N, nc);
%             dhi = zeros(N, nc);
%             while ~converged && min(alpha) > 0.01
%                 itCounter = 0;
%                 while itCounter < itLimit && ~converged
%                     % Finds hi for each component
%                     for i=1:nc
%                         hi(:,i) = (z(:,i) .* k(:,i)) ./ (fv .* (k(:,i) - 1) + 1);
%                         % Finds the derivative of hi for each component
%                         dhi(:,i) = (z(:,i) .* (k(:,i) - 1).^2) ./ ((fv .* (k(:,i) - 1) + 1).^2);
%                     end
%                     h = sum(hi, 2) - 1;
%                     dh = - sum(dhi, 2);
%                     
%                     % Update fv
%                     h(TwoPhase == 0) = 0;
%                     dh(TwoPhase == 0) = 1;
%                     fvnew = alpha .* (-h ./ dh) + fv;
%                     
%                     fv = fvnew;
%                     if norm(h, inf) < 1e-10
%                         converged = 1;
%                         disp(['Rachford-Rice converged in ', num2str(itCounter + 1), ' iterations, with alpha ', num2str(min(alpha))])
%                     end
%                     itCounter = itCounter + 1;
%                 end
%                 alpha (abs(h) > 1e-10) = alpha (abs(h) > 1e-10)/2;
%                 fv (fv > 1) = 0.9;
%                 fv (fv < 0) = 0.1;
%                 fv (isnan(fv)) = Status.Properties('ni_1').Value(isnan(fv));
%             end
%             if ~converged
%                 [~, cellIndex] = max(abs(h));
%                 disp('Warning: Flash did not fully converge!');
%                 disp(['The residual norm of the equilibrium equation is ', num2str(norm(h, inf)), ' in cell ', num2str(cellIndex)]);
%                 fv (fv > 1) = 1;
%                 fv(fv < 0 ) = 0;
%             end
            
            % find fv with bisection method
            converged = 0;
            itCounter = 0;
            hia = zeros(N, nc);
            hib = zeros(N, nc);
            hinew = zeros(N, nc);
            fva = 0.0 * ones(N, 1);
            fvb = 1 * ones(N, 1);
            while itCounter < itLimit && ~converged
                fvn = (fva + fvb)./2;
                % Finds hi for each component
                for i=1:nc
                    hia(:,i) = (z(:,i) .* k(:,i)) ./ (fva .* (k(:,i) - 1) + 1);
                    hib(:,i) = (z(:,i) .* k(:,i)) ./ (fvb .* (k(:,i) - 1) + 1);
                    hinew(:, i) = (z(:,i) .* k(:,i)) ./ (fvn .* (k(:,i) - 1) + 1);
                end
                ha = sum(hia, 2) - 1;
                hb = sum(hib, 2) - 1;
                hnew = sum(hinew, 2) - 1;
                % Update fv
                fva((ha .* hnew) > 0) = fvn((ha .* hnew) > 0);
                fvb((ha .* hnew) < 0) = fvn((ha .* hnew) < 0);
                hnew(TwoPhase == 0) = 0;
                if norm(hnew, inf) < 1e-10
                    converged = 1;
                    disp(['Rachford-Rice converged in ', num2str(itCounter + 1), ' iterations, with alpha ', num2str(min(alpha))])
                end
                itCounter = itCounter + 1;
            end
            fvn(SinglePhase == 1) = 1;
            fvn(SinglePhase == 2) = 0;
            disp(max(abs(fvn - fv)));
            
            %% 5. Solve for xs and ys
            % Have to make it general for nc components. Should be easy.
%             for i=1:nc
%                 x(TwoPhase==1, )
%             end
            % Solves for mole fractions in liquid phase
            x(TwoPhase == 1, 2) = z(TwoPhase == 1, 1) ./ (fv(TwoPhase == 1, 1) .* (k(TwoPhase == 1, 1) - 1) + 1);
            % Solves for mole fractions in gas phase
            x(TwoPhase == 1, 1) = k(TwoPhase == 1, 1) .* x(TwoPhase == 1, 2);
            x(:, 3:4) = 1 - x(:,1:2);
            
            Ni(:,1) = fv;
            Ni(:,2) = 1-fv;
            % Copy it into Status object
            for j=1:nph
                NI = Status.Properties(['ni_', num2str(j)]);
                NI.Value = Ni(:,j);
                for i=1:nc
                    X = Status.Properties(['x_', num2str(i),'ph',num2str(j)]);
                    X.Value = x(:,(i-1)*nph + j);
                end
            end
        end
    end
end