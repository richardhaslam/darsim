% Pressure convergence checker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_pressure < convergence_checker
    properties
        adm = 0
        OperatorsAssembler
        Incompressible = 0
    end
    methods 
        function PrintTitles(obj)
            disp(['Initial residual norm: ', num2str(obj.FirstResidualNorm, '%5.5e')]);
            disp('');
            disp('        ||Residual||   ||delta p||');
        end
        function converged = Check(obj, iter, residual, delta, Formulation, DiscretizationModel, State, LinearSolver)
            % Initialize
            converged = 0;
            % Compute Norms
            Norm2 = norm(delta, inf);
            if obj.adm
               R = LinearSolver.R;
               residual_c = R * residual;
               residual_c = residual_c ./ sum(R, 2);
               Norm1 = norm(residual_c, inf);
            else
               Norm1 = norm(residual, inf);
            end
            
            disp(['Iter ' num2str(iter, '%02d') '    ' num2str(Norm1, '%5.5e'), '    ', num2str(Norm2,'%5.5e')]);
            
            %Check convergence
            if (Norm1 < obj.ResidualTol || obj.Incompressible == 1) 
                converged = 1;
            end
        end
    end
end