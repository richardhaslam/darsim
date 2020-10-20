% Norm calculator geothermal_singlephase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 October 2016
%Last modified: 13 October 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef norm_calculator_geothermal_singlephase < norm_calculator
    properties
    end
    methods
        function [ResidualNorm, RHSNorm] = CalculateResidualNorm(obj, residual, RHS, N, Formulation)
            ResidualNorm = zeros(length(obj.FirstResidualNorm),1);
            ResidualNorm(1) = norm(residual(1 : N      ), 2);
            ResidualNorm(2) = norm(residual(N+1 : 2*N  ), 2);
            
            RHSNorm = zeros(length(obj.FirstRHSNorm),1);
            RHSNorm(1) = norm(RHS(1 : N      ), 2);
            RHSNorm(2) = norm(RHS(N+1 : 2*N  ), 2);
        end
        function [dp, dT] = CalculateSolutionNorm(obj, delta, N, State)
            dp = norm(delta(1:N)    ,2)/max(State.Properties('P_1').Value);
            dT = norm(delta(N+1:2*N),2)/max(State.Properties('T'  ).Value);
        end
    end
end