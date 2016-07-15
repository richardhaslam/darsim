% TimeLoop driver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef TimeLoop_Driver < handle
    properties
        TotalTime
        TStops %Choose with what frequency the solution as to be outputted
        MaxNumberOfTimeSteps
        Ndt
        dt
        Coupling
    end
    methods
        function obj = TimeLoop_Driver(total_time, n_reports)
            obj.TotalTime = total_time;
            obj.TStops = linspace(obj.TotalTime/n_reports, T, n_reports);
        end
        function AddCouplingStrategy(obj, coupling)
            obj.Coupling = coupling;
        end
        function Summary = SolveTimeDependentProblem(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, Summary, Writer)
            %%%%% START THE TIME LOOP %%%%%
            index = 1;         
            while (obj.TotalTime < obj.FinalTime && obj.obj.Ndt <= obj.MaxNumberOfTimeSteps)
                %% Initialise time-step
                disp(['Time-step ' num2str(obj.Ndt) ': Initial time: ' num2str(t/(3600*24),4) ' days']);
                tstart = tic;
                maxdT = obj.Tstops(index) - t;
                
                %% Non-linear Solver
                [ProductionSystem, obj.dt] =...
                    obj.Coupling.SolveTimeStep(ProductionSystem, FluidModel, DiscretizationModel, Formulation, maxdT);
                
                % Save Stats
                Summary = obj.Coupling.UpdateSummary(Summary, obj.Ndt);
                
                % Check for convergence at the end of the timestep
                if (obj.Coupling.NLSolver.Converged==0)
                   disp(['The solution has not converged at timestep ' num2str(obj.Ndt)]);
                   break
                end
                
                %% %%%Increase time and timestep counter
                disp('-----------------------------------------------')
                disp(['Final time: ' num2str((t+dT)/(3600*24),4) ' days, dT= ' num2str(dT) ' s']);
                disp(['end of time-step ' num2str(obj.Ndt)]);
                disp(char(5));
                t = t+dT;
                obj.Ndt = obj.Ndt + 1;
                
                % Saves the total timer of the timestep in the run summary
                Summary.Time(obj.Ndt-1) = t/(24*3600);
                Summary.NumberTimeSteps = obj.Ndt - 1;
                Summary.CouplingStats.SaveTimeStepTimer(obj.Ndt - 1, toc(tstart));
                
                %%%%%%%%%%%%%%PLOT SOLUTION%%%%%%%%%%%%%
                if (t == Tstops(index))
                    disp(['Printing solution to file at  ' num2str((t)/(3600*24),4) ' days'])
                    Writer.Plotter.PlotSolution(Status, Grid);
                    index = index + 1;
                end
            end
        end
    end
end