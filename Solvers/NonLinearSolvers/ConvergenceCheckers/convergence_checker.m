% Convergence checker for NLSolver base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker < handle
    properties
        Tol
        NormCalculator
        OldResidual
    end
    methods (Abstract)
        obj = Check(obj)
    end
    methods
        function output = Stagnating(obj, residual)
            output = 0;
            diff = abs(residual - obj.OldResidual)/abs(residual);
            obj.OldResidual = residual;
            if  diff < 1e-2
                output = 0;
            end
        end
    end
end