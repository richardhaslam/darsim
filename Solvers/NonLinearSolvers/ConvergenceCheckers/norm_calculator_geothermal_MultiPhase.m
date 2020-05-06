% Norm calculator comp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 October 2016
%Last modified: 13 October 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef norm_calculator_geothermal_MultiPhase < norm_calculator
    properties
        
    end
    methods
        function ResidualNorm = CalculateResidualNorm(obj, residual, N, Formulation)
            ResidualNorm = zeros(length(obj.FirstResidualNorm),1);
            ResidualNorm(1) = norm(residual(1 : N      ), 2);
            ResidualNorm(2) = norm(residual(N+1 : 2*N  ), 2);
        end
        function [dp, dT] = CalculateSolutionNorm(obj, delta, N, State)
            dp = norm(delta(1:N))/max(State.Properties('P_1').Value);
            dT = norm(delta(N+1:2*N))/max(State.Properties('hTfluid').Value);
        end
    end
end