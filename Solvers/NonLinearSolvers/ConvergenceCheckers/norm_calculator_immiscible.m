% Norm calculator immiscible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef norm_calculator_immiscible < norm_calculator
    methods
        function [ResidualNorm, RHSNorm, Equilibrium] = CalculateResidualNorm(obj, residual, RHS, N, Formulation)
            ResidualNorm = zeros(length(obj.FirstResidualNorm),1);
            ResidualNorm(1) = norm(residual(1 : N      ), inf);
            ResidualNorm(2) = norm(residual(N+1 : 2*N  ), inf);
            RHSNorm = zeros(length(obj.FirstRHSNorm),1);
            RHSNorm(1) = norm(RHS(1 : N      ), inf);
            RHSNorm(2) = norm(RHS(N+1 : 2*N  ), inf);
            Equilibrium = 0;
        end
        function [dp, dS] = CalculateSolutionNorm(obj, delta, N, State)
            dp = norm(delta(1:N), 2)/max(State.Properties('P_1').Value);
            dS = norm(delta(N+1:2*N), 2);
            % RelChange = delta(N+1:2*N) ./ State.Properties('S_1').Value;
            % RelChange(State.Properties('S_1').Value < 1e-4) = 0;
        end
    end
end