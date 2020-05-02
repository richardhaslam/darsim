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
        FirstResidual
        OldResidual
        FirstResidualNorm
        ResidualNorm
        FirstRHS
        FirstRHSNorm
        RHSNorm
    end
    methods (Abstract)
        obj = Check(obj) 
    end
    methods
        function ComputeFirstResidualNorm(obj, Residual, RHS, DiscretizationModel, LinearSolver)
            Nt = DiscretizationModel.N; 
            obj.FirstResidual = Residual;
            obj.FirstRHS      = RHS;
            % Compute Norms
            obj.FirstResidualNorm = zeros(obj.NumberOfEq,1);
            obj.FirstRHSNorm = zeros(obj.NumberOfEq,1);
            for eq = 1 : obj.NumberOfEq
                obj.FirstResidualNorm(eq) = norm(obj.FirstResidual((eq -1)*Nt+1:eq*Nt), inf); % The inf norm is better for LTS.
                obj.FirstRHSNorm(eq)      = norm(obj.FirstRHS(     (eq -1)*Nt+1:eq*Nt), inf);
            end
            obj.NormCalculator.FirstResidualNorm = obj.FirstResidualNorm;
            obj.NormCalculator.FirstRHSNorm      = obj.FirstResidualNorm;
            % Initializing the "ResidualNorm" as empty, to use it clean for
            % the coming iteration loop.
            obj.ResidualNorm = [];
            obj.RHSNorm = [];
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