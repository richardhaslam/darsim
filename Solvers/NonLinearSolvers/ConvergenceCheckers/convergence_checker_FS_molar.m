% Convergence checker for FS NLSolver molar formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_FS_molar < convergence_checker_FS
    properties
    end
    methods 
        function PrintTitles(obj, Residual)
            disp(['Initial residual norm: ', num2str(norm(Residual, inf))]);
            disp('');
            disp('        ||Residual||   ||Equilibrium||   ||delta p||   ||delta z||');
        end
        
    end
end