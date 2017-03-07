% Rachford-Rice Flash Calculator base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 12 January 2017
%Last modified: 12 January 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Standard_flash_calculator < Kvalues_flash_calculator
    properties
        tol = 1e-7;
        itMax = 1000;
    end
    methods
        function SinglePhase = Flash(obj, Status, Components, Phases)
            % Define SinglePhase objects
            nc = length(Components);
            nph = length(Phases);
            N = length(Status.Properties('z_1').Value);
            z = zeros(N, nc);
            Ni = zeros(N, nph);
            x = zeros(N, nc*nph);
            SinglePhase = zeros(N, 1);
            
            % Copy values in local variables
            for i=1:nc
                z(:, i) = Status.Properties(['z_', num2str(i)]).Value;
                for j=1:nph
                    x(:,(i-1)*nph + j) = Status.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
                end
            end
            for j=1:nph
                Ni(:,i) = Status.Properties(['ni_', num2str(i)]).Value;
            end
            ni = Ni(:,1);
            
            %% 1. Compute K-values
            k = obj.KvaluesCalculator.Compute(Status, Components, Phases);
            
            %% 2 Chek if we are in 2 phase region
            % 2.a: checking if it 's all liquid: checks if mix is below bubble
            BubCheck = z .* k;
            BubCheck = sum(BubCheck, 2);
            SinglePhase(BubCheck - 1 <= 0) = 2; % It s all liquid
            
            % 2.b: checking if it 's all vapor: checks if mix is above dew
            % point
            DewCheck = z ./ k;
            DewCheck = sum(DewCheck, 2);
            SinglePhase(DewCheck - 1 < 0) = 1;  % It s all vapour
            
            % Solve with Newton
            for i=1:N
                switch SinglePhase(i)
                    case(0)
                        converged = 0;
                        itcount = 1;
                        while ~converged && itcount < obj.itMax
                            Residual = zeros(2*nc+1,1);
                            J = zeros(2*nc+1);
                            for c=1:nc
                                Residual(c) = x(i, (c-1)*2+1) - k(i, c) * x(i, (c-1)*2+2);
                                Residual(nc+c) = z(i, c) - ni(i)*x(i, (c-1)*2+1) -(1-ni(i))*x(i, (c-1)*2+2);
                                Residual(2*nc+1) = Residual(2*nc+1) + x(i, (c-1)*2+1) - x(i, (c-1)*2+2);
                                % x_cv - k_c x_cl = 0
                                J (c, (c-1)*2+1:(c-1)*2+2) = [1, -k(i,c)];
                                % z_c - ni*x_cv - (1-ni)*x_cl = 0
                                J (nc + c, (c-1)*2+1:(c-1)*2+2) = [-ni(i), ni(i)-1];
                                J (nc + c, end) = x(i, (c-1)*2+2) - x(i, (c-1)*2+1);
                                % Sum x_cv-x_cl = 0
                                J (2*nc + 1, (c-1)*2+1:(c-1)*2+2) = [1, -1];
                            end
                            dx = -J\Residual;
                            x(i, :) = x(i, :) + dx(1:end-1)';
                            ni(i) = ni(i) + dx(2*nc+1);
                            if norm(Residual, inf) < obj.tol
                                converged = 1;
                            end
                            itcount = itcount+1;
                        end
                        if ~converged
                            disp(['Warning: Flash did not fully converge! in cell ', num2str(i)]);
                        end
                    case(1)
                        % It s all vapour
                        for c=1:nc
                            x(i,(c-1)*2+1) = z(i, c);
                            x(i,(c-1)*2+2) = 1; % This is for Stability
                        end
                        ni(i) = 1;
                    case(2)
                        % It s all liquid
                        for c = 1:nc
                            x(i, (c-1)*2+2) = z(i,c);
                            x(i, (c-1)*2+1) = 1; % This is for Stability
                        end
                        ni(i) = 0;
                end
            end
            Ni(:,1) = ni;
            Ni(:,2) = 1-ni;
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