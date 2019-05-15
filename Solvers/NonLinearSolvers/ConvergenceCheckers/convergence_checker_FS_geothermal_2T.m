% Convergence checker for Thermal simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_FS_geothermal_2T < convergence_checker
    properties
        AveragedTemperature
    end
    methods
        function ComputeFirstResidualNorm(obj, Residual, DiscretizationModel, LinearSolver)
            Nt = DiscretizationModel.N; 
            obj.FirstResidual = Residual;
            if obj.AveragedTemperature == "On"
                % Temperature Tf and Tr will be averaged and two equations will be summed
                Nm = DiscretizationModel.ReservoirGrid.N;
                obj.FirstResidual = zeros( (obj.NumberOfEq-1)*Nt+Nm , 1 );
                obj.FirstResidual(1:Nt) = Residual(1:Nt);
                obj.FirstResidual(Nt+1:Nt+Nm) = Residual(Nt+1:Nt+Nm)/2;
                obj.FirstResidual(Nt+Nm+1:2*Nt) = Residual(Nt+Nm+1:2*Nt);
                obj.FirstResidual(2*Nt+1:end) = Residual(Nt+1:Nt+Nm)/2;
            end
            % Compute Norms
            obj.FirstResidualNorm = zeros(obj.NumberOfEq,1);
            for eq = 1 : obj.NumberOfEq-1
                obj.FirstResidualNorm(eq) = norm(obj.FirstResidual((eq-1)*Nt+1:eq*Nt), 2);
            end
            obj.FirstResidualNorm(end) = norm(obj.FirstResidual(eq*Nt+1:end), 2);
            obj.NormCalculator.FirstResidualNorm = obj.FirstResidualNorm;
        end
        function PrintTitles(obj)
            disp(['Initial Fluid Pressure    Residual Norm: ', num2str(obj.FirstResidualNorm(1), '%5.5e')]);
            disp(['Initial Fluid Temperature Residual Norm: ', num2str(obj.FirstResidualNorm(2), '%5.5e')]);
            disp(['Initial Rock  Temperature Residual Norm: ', num2str(obj.FirstResidualNorm(3), '%5.5e')]);
            disp('');
            disp('           ||Residual P||   ||Residual Tf||   ||Residual Tr||   ||delta P||   ||delta Tf||   ||delta Tr|| ');
        end
        function converged = Check(obj, iter, residual, delta, Formulation, DiscretizationModel, State, LinearSolver)
            % Initialize
            converged = 0;
            Nt = DiscretizationModel.N;
            if obj.AveragedTemperature == "On"
                % Temperature Tf and Tr will be averaged and two equations will be summed
                Nm = DiscretizationModel.ReservoirGrid.N;
                residual_avg = residual;
                residual = zeros( (obj.NumberOfEq-1)*Nt+Nm , 1 );
                residual(1:Nt) = residual_avg(1:Nt);
                residual(Nt+1:Nt+Nm) = residual_avg(Nt+1:Nt+Nm)/2;
                residual(Nt+Nm+1:2*Nt) = residual_avg(Nt+Nm+1:2*Nt);
                residual(2*Nt+1:end) = residual_avg(Nt+1:Nt+Nm)/2;
                
                delta_ave = delta;
                delta = zeros( (obj.NumberOfEq-1)*Nt+Nm , 1 );
                delta(1:Nt) = delta_ave(1:Nt);
                delta(Nt+1:2*Nt) = delta_ave(Nt+1:2*Nt);
                delta(2*Nt+1:2*Nt+Nm) = delta_ave(Nt+1:Nt+Nm);
            end
            
            % Compute Norms
            [ResidualNorm] = obj.NormCalculator.CalculateResidualNorm(residual, Nt, Formulation);
            [dp, dTf, dTr] = obj.NormCalculator.CalculateSolutionNorm(delta, DiscretizationModel.N, State);
            disp(['Iter ', num2str(iter, '%02d') '-->   ', num2str(ResidualNorm(1), '%5.5e'), '      ' ...
                                                         , num2str(ResidualNorm(2), '%5.5e'), '      ' ...
                                                         , num2str(ResidualNorm(3), '%5.5e'), '      ' ...
                                                         , num2str(dp, '%5.5e'), '    ', num2str(dTf, '%5.5e'), '    ', num2str(dTr, '%5.5e')]);
            
            % check if it is stagnating
            stagnating = obj.Stagnating(ResidualNorm./obj.FirstResidualNorm);
            
            % Check convergence
            if ( (ResidualNorm(1) < obj.ResidualTol(1)) || (ResidualNorm(1)/obj.FirstResidualNorm(1) < obj.ResidualTol(1)) ) && ...
               ( (ResidualNorm(2) < obj.ResidualTol(2)) || (ResidualNorm(2)/obj.FirstResidualNorm(2) < obj.ResidualTol(2)) ) && ...
               ( (ResidualNorm(3) < obj.ResidualTol(3)) || (ResidualNorm(3)/obj.FirstResidualNorm(3) < obj.ResidualTol(3)) ) && ...
               ( dp < obj.SolutionTol(1) && dTf < obj.SolutionTol(2) && dTr < obj.SolutionTol(3) )
                converged = 1;
            elseif( isnan(ResidualNorm(1)) || isnan(ResidualNorm(2)) || isnan(ResidualNorm(3)) || stagnating || isnan(dp) || isnan(dTf) || isnan(dTr) || ~isreal([dp;dTf;dTr]) )
                converged = -1;
            end
        end
    end
end