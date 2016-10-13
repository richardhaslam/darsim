% Norm calculator comp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 October 2016
%Last modified: 13 October 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef norm_calculator_comp < norm_calculator
    methods
        function [Balance, Equilibrium] = ResidualNorm(obj, residual, N, Formulation)
            n_comp = Formulation.NofComponents;
            Balance = norm(residual(1:N*n_comp), inf);
            Equilibrium = norm(residual(N*n_comp+1:end), inf);
        end
        function [dp, dS] = SolutionNorm(obj, delta, N, State)
            dp = norm(delta(1:N), inf)/max(State.p);
            dS = norm(delta(N+1:end), inf);
        end
    end
end