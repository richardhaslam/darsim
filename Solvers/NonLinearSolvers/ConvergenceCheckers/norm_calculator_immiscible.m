% Norm calculator immiscible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef norm_calculator_immiscible < norm_calculator
    methods
        function [ResidualNorm, Equilibrium] = CalculateResidualNorm(obj, residual, N, Formulation)
            ResidualNorm = zeros(length(obj.FirstResidualNorm),1);
            ResidualNorm(1) = norm(residual(1 : N      ), 2);
            ResidualNorm(2) = norm(residual(N+1 : 2*N  ), 2);
            Equilibrium = 0;
        end
        function [dp, dS] = CalculateSolutionNorm(obj, delta, N, State)
            dp = norm(delta(1:N), 2)/max(State.Properties('P_1').Value);
            dS = norm(delta(N+1:2*N), 2)/max(State.Properties('S_1').Value);
            % RelChange = delta(N+1:2*N) ./ State.Properties('S_1').Value;
            % RelChange(State.Properties('S_1').Value < 1e-4) = 0;
        end
    end
end