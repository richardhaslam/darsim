% TimeLoop driver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef TimeLoop_Driver < handle
    properties
        Time
        TStops % List of solution report times
        TotalTime
        Ndt
        dt =0;
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
				% Adding the following three properties to the
                % "TimeStepSelector" class and assigning them for smarter
                % timestep selection.
                obj.Coupling.TimeStepSelector.FirstReportDt = obj.TStops(1);
                obj.Coupling.TimeStepSelector.BeforePreviousDt = obj.Coupling.TimeStepSelector.PreviousDt;
                obj.Coupling.TimeStepSelector.PreviousDt = obj.dt;
                obj.Coupling.TimeStepSelector.Index = index;
                
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
				DT = obj.Coupling.TimeStepSelector.Sec2DHMS(obj.dt);
                disp('-----------------------------------------------')
                disp(['Final time: ' num2str((obj.Time)/(3600*24),4) ' days, dt= ' num2str(obj.dt) ' sec (', ...
                      num2str(DT.Days), ' days : ', num2str(DT.Hours), ' hrs : ', num2str(DT.Minutes), ' mins : ', num2str(DT.Seconds), ' sec)']);																																  
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
                    Writer.WriteSolutionOnFile(ProductionSystem);
                    Writer.Index = Writer.Index + 1;
                    index = index + 1;
                end
            end
        end
		function Time = Sec2DHMS(obj, TimeInSec)
            Time.Days  = floor( TimeInSec / (24*3600) );
            TimeInSec  = TimeInSec - Time.Days * 24*3600;
            
            Time.Hours = floor( TimeInSec / 3600 );
            TimeInSec  = TimeInSec - Time.Hours * 3600;
            
            Time.Minutes = floor( TimeInSec / 60 );
            
            Time.Seconds  = TimeInSec - Time.Minutes * 60;
        end
    end
end

