% Convergence checker for ADM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 2 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_ADM_geothermal_2T < convergence_checker_FS_geothermal_2T
    properties
        OperatorsAssembler
    end
    methods
        function ComputeFirstResidualNorm(obj, Residual, DiscretizationModel, LinearSolver)
            Nt_ADM = DiscretizationModel.ADMGrid.Ntot;
            % Get ADM Operators
            ADMRest = LinearSolver.R;
            % Restrict first residual
            Residual_ADM = ADMRest * Residual ./ sum(ADMRest, 2);
            obj.FirstResidual = Residual_ADM;
            if obj.AveragedTemperature == "On"
                % Temperature Tf and Tr will be averaged and two equations will be summed
                Nm_ADM = DiscretizationModel.ADMGrid_Reservoir.Ntot;
                obj.FirstResidual = zeros( (obj.NumberOfEq-1)*Nt_ADM+Nm_ADM , 1 );
                obj.FirstResidual(1:Nt_ADM) = Residual_ADM(1:Nt_ADM);
                obj.FirstResidual(Nt_ADM+1:Nt_ADM+Nm_ADM) = Residual_ADM(Nt_ADM+1:Nt_ADM+Nm_ADM)/2;
                obj.FirstResidual(Nt_ADM+Nm_ADM+1:2*Nt_ADM) = Residual_ADM(Nt_ADM+Nm_ADM+1:2*Nt_ADM);
                obj.FirstResidual(2*Nt_ADM+1:end) = Residual_ADM(Nt_ADM+1:Nt_ADM+Nm_ADM)/2;
            end
            % Compute Norms
            obj.FirstResidualNorm = zeros(obj.NumberOfEq,1);
            for eq = 1 : obj.NumberOfEq-1
                obj.FirstResidualNorm(eq) = norm(obj.FirstResidual((eq-1)*Nt_ADM+1:eq*Nt_ADM), 2);
            end
            obj.FirstResidualNorm(end) = norm(obj.FirstResidual(eq*Nt_ADM+1:end), 2);
            obj.NormCalculator.FirstResidualNorm = obj.FirstResidualNorm;
        end
        function converged = Check(obj, iter, residual, delta, Formulation, DiscretizationModel, State, LinearSolver)
            % Initialize
            converged = 0;
            Nt = DiscretizationModel.N;
            Nt_ADM = DiscretizationModel.ADMGrid.Ntot;
            
            ADMRest = LinearSolver.R;  % Get ADM Operators
            Residual_ADM = ADMRest * residual ./ sum(ADMRest, 2); % Restrict each residual and divide by number of cells in each coarse node
            
            if obj.AveragedTemperature == "On"
                Nm = DiscretizationModel.ReservoirGrid.N;
                Nm_ADM = DiscretizationModel.ADMGrid_Reservoir.Ntot;
                Residual_ADM_Avg = Residual_ADM;
                Residual_ADM = zeros( (obj.NumberOfEq-1)*Nt_ADM+Nm_ADM , 1 );
                Residual_ADM(1:Nt_ADM) = Residual_ADM_Avg(1:Nt_ADM);
                Residual_ADM(Nt_ADM+1:Nt_ADM+Nm_ADM) = Residual_ADM_Avg(Nt_ADM+1:Nt_ADM+Nm_ADM)/2;
                Residual_ADM(Nt_ADM+Nm_ADM+1:2*Nt_ADM) = Residual_ADM_Avg(Nt_ADM+Nm_ADM+1:2*Nt_ADM);
                Residual_ADM(2*Nt_ADM+1:end) = Residual_ADM_Avg(Nt_ADM+1:Nt_ADM+Nm_ADM)/2;
                
                delta_ave = delta;
                delta = zeros( (obj.NumberOfEq-1)*Nt+Nm , 1 );
                delta(1:Nt) = delta_ave(1:Nt);
                delta(Nt+1:2*Nt) = delta_ave(Nt+1:2*Nt);
                delta(2*Nt+1:2*Nt+Nm) = delta_ave(Nt+1:Nt+Nm);
            end

            % Compute Norms            
            [ResidualNorm] =  obj.NormCalculator.CalculateResidualNorm(Residual_ADM, Nt_ADM, Formulation);
            [dp, dTf, dTr] = obj.NormCalculator.CalculateSolutionNorm(delta, Nt, State);
            
            disp(['Iter ', num2str(iter, '%02d') '-->   ', num2str(ResidualNorm(1)/obj.FirstResidualNorm(1), '%5.5e'), '      ' ...
                                                         , num2str(ResidualNorm(2)/obj.FirstResidualNorm(2), '%5.5e'), '      ' ...
                                                         , num2str(ResidualNorm(3)/obj.FirstResidualNorm(3), '%5.5e'), '      ' ...
                                                         , num2str(dp, '%5.5e'), '    ', num2str(dTf, '%5.5e'), '    ', num2str(dTr, '%5.5e')]);
            
            % check if it is stagnating
            stagnating = obj.Stagnating(ResidualNorm./obj.FirstResidualNorm);
            
            % Check convergence
            if ( (ResidualNorm(1) < obj.ResidualTol(1)) || (ResidualNorm(1)/obj.FirstResidualNorm(1) < obj.ResidualTol(1)) ) && ...
               ( (ResidualNorm(2) < obj.ResidualTol(2)) || (ResidualNorm(2)/obj.FirstResidualNorm(2) < obj.ResidualTol(2)) ) && ...
               ( (ResidualNorm(3) < obj.ResidualTol(3)) || (ResidualNorm(3)/obj.FirstResidualNorm(3) < obj.ResidualTol(3)) ) && ...
               ( dp < obj.SolutionTol(1) && dTf < obj.SolutionTol(2) && dTr < obj.SolutionTol(3) )
                converged = 1;
            elseif(isnan(ResidualNorm(1)) || isnan(ResidualNorm(2)) || isnan(ResidualNorm(3)) || stagnating || isnan(dp) || isnan(dTf) || isnan(dTr))
                converged = -1;
            end
        end
    end
end