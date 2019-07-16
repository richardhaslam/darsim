% Convergence checker for NLSolver base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker < handle
    properties
        NumberOfEq
        ResidualTol
        SolutionTol
        NormCalculator
        FirstResidualNorm
        FirstResidual
        OldResidual
    end
    methods (Abstract)
        obj = Check(obj) 
    end
    methods
        function ComputeFirstResidualNorm(obj, Residual, DiscretizationModel, LinearSolver)
            Nt = DiscretizationModel.N; 
            obj.FirstResidual = Residual;
            % Compute Norms
            obj.FirstResidualNorm = zeros(obj.NumberOfEq,1);
            for eq = 1 : obj.NumberOfEq
                obj.FirstResidualNorm(eq) = norm(obj.FirstResidual((eq -1)*Nt+1:eq*Nt), inf); % The inf norm is better for LTS.
            end
            obj.NormCalculator.FirstResidualNorm = obj.FirstResidualNorm;
        end
        function output = Stagnating(obj, residual, delta)
            output = 0;
            if isempty(obj.OldResidual)
                obj.OldResidual = zeros(length(residual),1);
            end
            diff = abs(residual - obj.OldResidual)./abs(residual);
            obj.OldResidual = residual;
            if  diff < 1e-2
                output = 0;
            end
        end
    end
end