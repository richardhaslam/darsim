% TimeLoop driver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 18 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef TimeLoop_Driver < handle
    properties
        TotalTime
        Time
        TStops % List of solution report times
        MaxNumberOfTimeSteps
        Ndt
        dt
        Coupling
    end
    methods
        function obj = TimeLoop_Driver(total_time, n_reports)
            obj.TotalTime = total_time;
            obj.TStops = linspace(obj.TotalTime/n_reports, obj.TotalTime, n_reports);
            % obj.TStops = [1.9602e+06*0.051; 1.9602e+06*0.51; 1.9602e+06*2.55; 1.9602e+06*5.1; obj.TotalTime];
        end
        function AddCouplingStrategy(obj, coupling)
            obj.Coupling = coupling;
        end
        function Summary = SolveTimeDependentProblem(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, Summary, Writer)
            %%%%% START THE TIME LOOP %%%%%
            index = 1;   
            obj.Time = 0;
            obj.Ndt = 1;
            while (obj.Time < obj.TotalTime && obj.Ndt <= obj.MaxNumberOfTimeSteps)
                %% Initialise time-step
                disp(['Time-step ' num2str(obj.Ndt) ': Initial time: ' num2str(obj.Time/(3600*24),4) ' days']);
                tstart = tic;
                
                %% Non-linear Solver
                obj.Coupling.TimeStepSelector.ReportDt = obj.TStops(index) - obj.Time;

                obj.dt = obj.Coupling.SolveTimeStep(ProductionSystem, FluidModel, DiscretizationModel, Formulation);                
                
                % Average for ADM
                DiscretizationModel.AverageMassOnCoarseBlocks(ProductionSystem.Reservoir.State, FluidModel, Formulation);
                
                % Save Stats
                obj.Coupling.UpdateSummary(Summary, ProductionSystem.Wells, obj.Ndt, obj.dt);
                
                % Check for convergence at the end of the timestep
                if ~obj.Coupling.Converged
                   disp(['The solution has not converged at timestep ' num2str(obj.Ndt)]);
                   break
                end
                
                %% %%%Increase time and timestep counter
                obj.Time = obj.Time + obj.dt;
                disp('-----------------------------------------------')
                disp(['Final time: ' num2str((obj.Time)/(3600*24),4) ' days, dT= ' num2str(obj.dt) ' s']);
                disp(['end of time-step ' num2str(obj.Ndt)]);
                disp(char(5));
                obj.Ndt = obj.Ndt + 1;
                
                % Saves the total timer of the timestep in the run summary
                Summary.Time(obj.Ndt-1) = obj.Time/(24*3600);
                Summary.NumberTimeSteps = obj.Ndt - 1;
                Summary.CouplingStats.SaveTimeStepTimer(obj.Ndt - 1, toc(tstart));
                Summary.SaveGridStats(obj.Ndt - 1, DiscretizationModel);
                
                %%%%%%%%%%%%%%PLOT SOLUTION%%%%%%%%%%%%%
                if (obj.Time == obj.TStops(index))
                    disp(['Printing solution to file at  ' num2str((obj.Time)/(3600*24),4) ' days'])
                    disp(char(5));
                    Writer.PlotSolution(ProductionSystem, DiscretizationModel);
                    index = index + 1;
                end
            end
        end
    end
end