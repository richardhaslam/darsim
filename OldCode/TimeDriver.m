%        Time loop driver             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last modified: 15 May 2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time loop driver     
%The time-dependent problem is solved in this script.
%Two solution strategies exist:
%                       1. Sequential
%                           a.IMPES
%                           b.IMPSAT
%                       2. FIM
%Variables:
%   T: total simulation time
%   TimeSteps: maximum number of timesteps allowed
%   Grid: contains all fine-scale grid information
%   CoarseGrid: contains ADM coarse grids information
%   Fluid: fluids properties
%   FIM: settings of FIM non-linear solver
%   Sequential: settings of Sequential non-linear solver
%   ADMSettings: settings of ADM solver
%   Inj: injection wells
%   Prod: production wells

%   Status: contains the status of the reservoir (p, s, pc, x)

%   dT: timestep size
%   Ndt: timestep number
%   t: current simulation time


function [Status, Summary] = ...
    TimeDriver...
    (Grid, K, Fluid, Status, Inj, Prod, T, TimeSteps, Summary, Coupling, ...
    Output, CoarseGrid, ADMSettings, FlashSettings)
%%%%% START THE TIME LOOP %%%%%

t = 0;
Ndt = 1; 
Converged = 0;
index = 1;
%Choose with what frequency the solution as to be outputted
Tstops = linspace(T/10, T, 10);

while (t<T && Ndt <= TimeSteps)
    %% Initialise time-step
    disp(['Time-step ' num2str(Ndt) ': Initial time: ' num2str(t/(3600*24),4) ' days']);
    tstart = tic;

    maxdT = Tstops - t;
    %Reset all fine cells to active
    Grid.Active = ones(Grid.N, 1);
    
    %% Non-linear Solver
    switch (Coupling.Name)
        case ('Sequential')
            disp('--------------Sequential non-linear solver-----');
            if ~Sequential.ImpSat
                Sequential.MaxExtIter = 1;
            end
            [Status, Inj, Prod, dT, Converged, Timers, Sequential.ImplicitSolver] =...
                SequentialStrategy(Status, K, Grid, Fluid, Inj, Prod, Sequential, Ndt, maxdT(index));
        case ('FIM')
            disp('------------FIM Non-linear solver--------------');
            if (Ndt==1)
                % Use IMPES to estimate timestep size for the 1st timestep
                [~, U, ~, ~, ~, ~] = PressureSolver(Grid, Inj, Prod, Fluid, Status.s, K);
                dTnext = timestepping(Fluid, Grid, Coupling.CFL, U);
                dT = min(dTnext, maxdT(index));
                
                [Status, dT, dTnext, Inj, Prod, Summary, Converged, CoarseGrid, Grid] = ...
                    FIMNonLinearSolver...
                (Status, K, Grid, Fluid, Inj, Prod, Coupling.NLSolver, Summary, dT, Ndt, CoarseGrid, ADMSettings, FlashSettings);
            else
                dT = min(dTnext, maxdT(index));

                %Non-linear solver
                [Status, dT, dTnext, Inj, Prod, Summary, Converged, CoarseGrid, Grid] = ...
                    FIMNonLinearSolver...
                (Status, K, Grid, Fluid, Inj, Prod, Coupling.NLSolver, Summary, dT, Ndt, CoarseGrid, ADMSettings, FlashSettings);
            end
    end
    
    %Check for convergence at the end of the timestep
    if (Converged==0) 
        switch(Strategy) 
            case('Sequential')
                disp(['The solution has not converged at timestep ' num2str(Ndt)]);
                break
            case('FIM')
                disp(['The solution has not converged at time-step '...
                    num2str(Ndt) ' even with several time-step chops']);
                break
        end
    end
    %%end of Non-linear Solver
    
    %% %%%Increase time and timestep counter
    disp('-----------------------------------------------')
    disp(['Final time: ' num2str((t+dT)/(3600*24),4) ' days, dT= ' num2str(dT) ' s']);
    disp(['end of time-step ' num2str(Ndt)]);
    disp(char(5));
    t = t+dT;    
    Ndt = Ndt + 1;
    
    % Saves the total timer of the timestep in the run summary
    Summary.Time(Ndt-1) = t/(24*3600);
    Summary.NumberTimeSteps = Ndt - 1;
    Summary.CouplingStats.SaveTimeStepTimer(Ndt - 1, toc(tstart)); 
    
    %%%%%%%%%%%%%%PLOT SOLUTION%%%%%%%%%%%%%
    if (t == Tstops(index))
        disp(['Printing solution to file at  ' num2str((t)/(3600*24),4) ' days'])
        Output.Plotter.PlotSolution(Status, Grid);
        index = index +1;
    end
end
end