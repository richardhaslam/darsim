% Norm calculator comp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef norm_calculator_compositional < norm_calculator
    methods
        function [Balance, RHSNorm, Equilibrium] = CalculateResidualNorm(obj, residual, RHS, N, Formulation)
            Balance = zeros(length(obj.FirstResidualNorm),1);
            RHSNorm = zeros(length(obj.FirstRHSNorm)     ,1);
            n_comp = Formulation.NofComponents;
            for i=1:n_comp
                Balance(i) = norm(residual((i-1)*N + 1:i*N), 2);
                RHSNorm(i) = norm(RHS     ((i-1)*N + 1:i*N), 2);
            end
            Equilibrium = norm(residual(N*n_comp+1:end), inf);
        end
        function [dp, dS] = CalculateSolutionNorm(obj, delta, N, State)
            dp = norm(delta(1:N), inf)/max(State.Properties('P_2').Value);
            dS = norm(delta(N+1:end), inf);
            %[a, b] = max(norm(dS));
            %disp(['max is ', num2str(a), ' in cell ', num2str(b)]);
        end
    end
end