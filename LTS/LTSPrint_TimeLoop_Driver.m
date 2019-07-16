% TimeLoop driver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTSPrint_TimeLoop_Driver <  TimeLoop_Driver

    methods
        function Summary = SolveTimeDependentProblem(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, Summary, Writer)
            %%%%% START THE TIME LOOP %%%%%
            index = 1;
            index_internal = 1;
            obj.Time = 0;
            obj.Ndt = 1;
            EndOfSimCriterion = 0;
            while EndOfSimCriterion == 0
                %% Initialise time-step
                disp(['Time-step ' num2str(obj.Ndt) ': Initial time: ' num2str(obj.Time/(3600*24),4) ' days']);
                tstart = tic;
                
                %% Solve Coupled problem at time-step n
                obj.Coupling.TimeStepSelector.ReportDt = obj.TStops(index) - obj.Time;
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
                disp('-----------------------------------------------')
                disp(['Final time: ' num2str((obj.Time)/(3600*24),4) ' days, dt= ' num2str(obj.dt) ' s']);
                disp(['end of time-step ' num2str(obj.Ndt)]);
                disp(newline);
                obj.Ndt = obj.Ndt + 1;

                % Saves the total time of the timestep in the run summary
                Summary.Time(obj.Ndt-1) = obj.Time/(24*3600);
                Summary.NumberTimeSteps = obj.Ndt - 1;
                Summary.CouplingStats.SaveTimeStepTimer(obj.Ndt - 1, toc(tstart));
                Summary.SaveGridStats(obj.Ndt - 1, DiscretizationModel);
%                Summary.SaveLTSiter( obj.Ndt, obj.Coupling.LTS_iters);
                
                %% Has simulation ended?
                EndOfSimCriterion = obj.EndOfSimEvaluator.HasSimulationEnded(EndOfSimCriterion, Summary, ProductionSystem, obj.Time, obj.Ndt);
                
                %% %%%%%%%%%%%%PLOT SOLUTION%%%%%%%%%%%%%
                  %%%%% PLOT LTS  NEW SOLUTION %%%%%%%%
                if (obj.Time == obj.TStops(index) || obj.TStops(index) == 0 ||EndOfSimCriterion==1 )
                    disp(['Printing solution to file at  ' num2str((obj.Time)/(3600*24),4) ' days'])
                    disp(newline);
                    if size(obj.Coupling.ActCellsSummary,2) ~= 0 && sum(sum(obj.Coupling.ActCellsSummary)) ~= 0
                        % for each sub/ref we modify the ActiveTime attribute
                        % and the ProductionSystem.State
                        % of Reservoir Grid to plot it.
                        StateGlob = status();
                        StateGlob.CopyProperties(ProductionSystem.Reservoir.State);
                        
                        for i = 1:size(obj.Coupling.ActCellsSummary,2)
                            DiscretizationModel.ReservoirGrid.ActiveTime = obj.Coupling.ActCellsSummary(:,i);
                            ProductionSystem.Reservoir.State.CopyProperties(obj.Coupling.StatesSummary(i));
                            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                            Writer.PlotSolution(ProductionSystem, DiscretizationModel);
                            Writer.WriteSolutionOnFile(ProductionSystem, index_internal);
                            index_internal = index_internal + 1;
                        end
                        DiscretizationModel.ReservoirGrid.ActiveTime = ones(size(obj.Coupling.ActCellsSummary(:,i)));
                        ProductionSystem.Reservoir.State.CopyProperties(StateGlob);
                        ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                        index = index + 1;
                    else
                    Writer.PlotSolution(ProductionSystem, DiscretizationModel);
                    Writer.WriteSolutionOnFile(ProductionSystem, index_internal);
                    index_internal = index_internal + 1;
                    index = index + 1;
                    end
                end
            end
        end
    end
end