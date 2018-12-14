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
        OldResidual
    end
    methods (Abstract)
        obj = Check(obj)
    end
    methods
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