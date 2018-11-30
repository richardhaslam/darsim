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
        Time
        TStops % List of solution report times
        TotalTime
        Ndt
        dt
        Coupling
        EndOfSimEvaluator
    end
    methods
        function obj = TimeLoop_Driver(n_reports, TotalTime, MaxNumTimesteps)
            obj.TotalTime = TotalTime;
            if n_reports == 0
                % Prnt at every time-step
                obj.TStops = TotalTime * zeros(MaxNumTimesteps ,1);
            else
                obj.TStops = linspace(TotalTime/n_reports, TotalTime, n_reports);
            end
            %obj.TStops = [365*24*3600:365*24*3600:1000*365*24*3600];
        end
        function AddCouplingStrategy(obj, coupling)
            obj.Coupling = coupling;
        end
        function AddEndOfSimEvaluator(obj, end_of_sim_eval)
            obj.EndOfSimEvaluator = end_of_sim_eval;
        end
        function Summary = SolveTimeDependentProblem(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, Summary, Writer)
            %%%%% START THE TIME LOOP %%%%%
            index = 1;   
            obj.Time = 0;
            obj.Ndt = 1;
            EndOfSimCriterion = 0;
            while EndOfSimCriterion == 0
                %% Initialise time-step
                disp(['Time-step ' num2str(obj.Ndt) ': Initial time: ' num2str(obj.Time/(3600*24),4) ' days']);
                tstart = tic;
                
                %% Solve Coupled problem at time-step n
                obj.Coupling.TimeStepSelector.ReportDt = obj.TStops(index) - obj.Time;

                [obj.dt, EndOfSimCriterion] = obj.Coupling.SolveTimeStep(ProductionSystem, FluidModel, DiscretizationModel, Formulation);                
                
                % Average for ADM
                DiscretizationModel.AverageMassOnCoarseBlocks(ProductionSystem, FluidModel, Formulation);
                
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
                disp(newline);
                obj.Ndt = obj.Ndt + 1;

                % Saves the total time of the timestep in the run summary
                Summary.Time(obj.Ndt-1) = obj.Time/(24*3600);
                Summary.NumberTimeSteps = obj.Ndt - 1;
                Summary.CouplingStats.SaveTimeStepTimer(obj.Ndt - 1, toc(tstart));
                Summary.SaveGridStats(obj.Ndt - 1, DiscretizationModel);
                
                %% Has simulation ended?
                EndOfSimCriterion = obj.EndOfSimEvaluator.HasSimulationEnded(EndOfSimCriterion, Summary, ProductionSystem, obj.Time, obj.Ndt);
                
                %% %%%%%%%%%%%%PLOT SOLUTION%%%%%%%%%%%%%
                if (obj.Time == obj.TStops(index) || obj.TStops(index) == 0 ||EndOfSimCriterion==1 )
                    disp(['Printing solution to file at  ' num2str((obj.Time)/(3600*24),4) ' days'])
                    disp(newline);
                    Writer.PlotSolution(ProductionSystem, DiscretizationModel);
                    Writer.WriteSolutionOnFile(ProductionSystem, index);
                    index = index + 1;
                end
            end
        end
    end
end