% Convergence checker for transport solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_transport < convergence_checker
    properties
        adm = 0
        lts = 0
        capillary = 0
        OperatorsAssembler
    end
    methods 
        function PrintTitles(obj)
            disp(['Initial residual norm: ', num2str(obj.FirstResidualNorm, '%5.5e')]);
            disp('');
            disp('        ||Residual||   ||delta S||');
        end
        function converged = Check(obj, iter, residual, delta, Formulation, DiscretizationModel, State, LinearSolver)
            % Initialize
            converged = 0;
            % Compute Norms
            if obj.adm  
                    R = LinearSolver.R;
                    residual_c = R * residual;
                    residual_c = residual_c ./ sum(R, 2);
                    Norm1 =  norm(residual_c, inf);
            else
                Norm1 =  norm(residual, inf);
            end
            Norm2 = norm(delta, inf);
            
            disp(['Iter ' num2str(iter) '    ' num2str(Norm1, '%5.5e'), '    ', num2str(Norm2,'%5.5e')]);
            
            %Check convergence
            if obj.capillary
                if (Norm1 < obj.ResidualTol)
                    converged = 1;
                end
            else
                if (Norm1 < obj.ResidualTol && Norm2 < obj.SolutionTol)
                converged = 1;
                end
            end
        end
        function converged = CheckLTS(obj, iter, residual, delta, Formulation, DiscretizationModel, State, LinearSolver, ActCells)
            % Initialize
            converged = 0;
            % Compute Norms
            if obj.adm
                if obj.lts
                    Rred = LinearSolver.Rred;
                    ActCells = logical(ActCells);
            
                    residual_c = Rred * residual(ActCells);
                    residual_c = residual_c ./ sum(Rred, 2);
                    Norm1 =  norm(residual_c, inf);
                else   
                    R = LinearSolver.R;
                    residual_c = R * residual;
                    residual_c = residual_c ./ sum(R, 2);
                    Norm1 =  norm(residual_c, inf);
                end
            else
                Norm1 =  norm(residual, inf);
            end
            Norm2 = norm(delta, inf);
            
            disp(['Iter ' num2str(iter) '    ' num2str(Norm1, '%5.5e'), '    ', num2str(Norm2,'%5.5e')]);
            
            %Check convergence
            if obj.capillary
                if (Norm1 < obj.ResidualTol)
                    converged = 1;
                end
            else
                if (Norm1 < obj.ResidualTol && Norm2 < obj.SolutionTol)
                    converged = 1;
                end
            end
        end
    end
end