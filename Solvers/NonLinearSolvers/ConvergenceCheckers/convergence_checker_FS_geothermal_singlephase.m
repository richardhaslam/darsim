% Convergence checker for Thermal simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_FS_geothermal_singlephase < convergence_checker
    properties
    end
    methods    
        function converged = CheckConvergenceCondition(obj,iter,dp,dT)
            disp(['Iter ', num2str(iter, '%02d') '-->   ', num2str(obj.ResidualNorm(iter,1), '%5.5e'), '      ' ...
                                                         , num2str(obj.ResidualNorm(iter,2), '%5.5e'), '      ' ...
                                                         , num2str(dp, '%5.5e'), '    ', num2str(dT, '%5.5e')]);
            % Initialize
            converged = 0;
            
            % check if it is stagnating
            stagnating = obj.Stagnating(obj.ResidualNorm(iter,:)./obj.FirstResidualNorm);
            
            %Check convergence
            if ( (obj.ResidualNorm(iter,1) < obj.ResidualTol(1)) || (obj.ResidualNorm(iter,1)/obj.FirstResidualNorm(1) < obj.ResidualTol(1)) ) && ...
               ( (obj.ResidualNorm(iter,2) < obj.ResidualTol(2)) || (obj.ResidualNorm(iter,2)/obj.FirstResidualNorm(2) < obj.ResidualTol(2)) ) && ...
               ( dp < obj.SolutionTol(1) && dT < obj.SolutionTol(2) )
                converged = 1;
            elseif iter>5 && dp < obj.SolutionTol(1) && dT < obj.SolutionTol(2) && ...
                   abs(obj.ResidualNorm(end-1,1) - obj.ResidualNorm(end-2,1)) < 1e-3 * obj.ResidualNorm(end-2,1) && ...
                   abs(obj.ResidualNorm(end  ,1) - obj.ResidualNorm(end-1,1)) < 1e-3 * obj.ResidualNorm(end  ,1) && ...
                   abs(obj.ResidualNorm(end-1,2) - obj.ResidualNorm(end-2,2)) < 1e-3 * obj.ResidualNorm(end-2,2) && ...
                   abs(obj.ResidualNorm(end  ,2) - obj.ResidualNorm(end-1,2)) < 1e-3 * obj.ResidualNorm(end  ,2)
               converged = 1;
            elseif (isnan(obj.ResidualNorm(iter,1)) || isnan(obj.ResidualNorm(iter,2)) || stagnating || isnan(dp) || isnan(dT))
                converged = -1;
            end
        end
    end
end