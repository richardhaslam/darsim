% Norm calculator immiscible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef norm_calculator_immiscible < norm_calculator
    methods
        function [Balance, Equilibrium] = ResidualNorm(obj, residual, N, Formulation)
            Balance = norm(residual, inf);
            Equilibrium = 0;
        end
        function [dp, dS] = SolutionNorm(obj, delta, N, State)
            dp = norm(delta(1:N), inf)/max(State.Properties('P_1').Value);
            % RelChange = delta(N+1:2*N) ./ State.Properties('S_1').Value;
            % RelChange(State.Properties('S_1').Value < 1e-4) = 0;
            dS = norm(delta(N+1:2*N), inf);
        end
    end
end