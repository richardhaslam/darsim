% Norm calculator immiscible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 October 2016
%Last modified: 13 October 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef norm_calculator_immiscible < norm_calculator
    methods
        function [Balance, Equilibrium] = ResidualNorm(obj, residual, N, Formulation)
            Balance = norm(residual, inf);
            Equilibrium = 0;
        end
        function [dp, dS] = SolutionNorm(obj, delta, N, State)
            dp = norm(delta(1:N), inf)/max(State.Properties('P_1').Value);
            dS = norm(delta(N+1:end), inf);
        end
    end
end