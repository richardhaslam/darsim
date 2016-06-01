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

%   P: pressure at current timestep
%   S: saturation at current timestep 
%   P0: pressure at previous timestep
%   S0: saturation at previous timestep
%   Pc: capillary pressure

%   dT: timestep size
%   Ndt: timestep number
%   t: current simulation time

%TIMERS:
%   TimerTimestep: total time of each timestep
%   For Sequential non-linear solver
%       TimerPressure: pressure solver timer
%       TimerBalance: balance check timer
%       TimerSaturation: saturation solver timer
%   For FIM non-linear solver
%       TimerConstruct: Jacobian construction timer
%       TimerSolve: Linear system solve timer
%   For ADM solver
%       TimerRP: construction of R and P operators

%% Timers and variables for statistics
TimerTimestep = zeros(TimeSteps,1);
if (strcmp(Strategy, 'Sequential')==1)
    TimerPressure = zeros(TimeSteps,1);
    TimerBalance = zeros(TimeSteps,1);
    TimerSaturation = zeros(TimeSteps,1);
else
    TimerConstruct = zeros(TimeSteps,1);
    TimerSolve = zeros(TimeSteps,1);
    TimerRP = zeros(TimeSteps,1);
end

%%%%% START THE TIME LOOP %%%%%

t = 0;
Ndt = 1; 
Converged = 0;
index = 1;
Saturations = zeros(Grid.N, 50);
Pressures = zeros(Grid.N, 50);
NwProduction = zeros(length(Prod) + 1, TimeSteps);
WProduction = zeros(length(Prod) + 1, TimeSteps);
vtkcount = 1;
%Choose with what frequency the solution as to be outputted
Tstops = linspace(T/50, T, 50);


while (t<T && Ndt <= TimeSteps)
    %% Initialise time-step
    disp(['Time-step ' num2str(Ndt) ': Initial time: ' num2str(t/(3600*24),4) ' days']);
    tstart = tic;
    S0 = S;
    P0 = P;
    maxdT = Tstops - t;
    %Reset all fine cells to active
    Grid.Active = ones(Grid.N, 1);
    
    %% Non-linear Solver
    switch (Strategy)
        case ('Sequential')
            disp('--------------Sequential non-linear solver-----');
            if ~Sequential.ImpSat
                Sequential.MaxExtIter = 1;
            end
            [P, S, Pc, Inj, Prod, dT, Converged, Timers, Sequential.ImplicitSolver] =...
                SequentialStrategy(S0, K, Grid, Fluid, Inj, Prod, Sequential, Ndt, maxdT(index));
        case ('FIM')
            disp('------------FIM Non-linear solver--------------');
            disp('        ||Residual||   Sat. delta');
            FIM.timestep (Ndt) = Ndt;
            if (Ndt==1)
                % Use IMPES as intial guess for pressure for the 1st timestep
                [P0, U, Pc, Wells, Inj, Prod] = PressureSolver(Grid, Inj, Prod, Fluid, S, K);
                %[Pms, ~] = MMsFVPressureSolver(Grid, Inj, Prod, K, Fluid, S, CoarseGrid, maxLevel);
                
                %Plot initial conditions
                switch (Options.PlotSolution)
                    case('Matlab')
                        if ADMSettings.active
                            Plotting_ADM
                        else
                            Plotting(Grid, P0, Pc, S, Fluid, 'red', 'blue', Prod, Inj);
                        end
                    case('VTK')
                        Wells2VTK(Grid, Inj, Prod, Directory, Problem);
                        Write2VTK(Directory, Problem, vtkcount, Grid, K, P0, S, Pc, ADMSettings.active, CoarseGrid, ADMSettings.maxLevel, 1);
                        vtkcount = vtkcount + 1;
                end
                
                %Keep first timestep to be small
                Grid.CFL = 1;
                maxiteration = FIM.MaxIter;
                FIM.MaxIter = 40;
                dT = timestepping(Fluid, S, Grid, U, Wells);

                %Compute rock transmissibility
                [Trx, Try] = ComputeTransmissibility(Grid, K);
                
                %Non-linear solver
                [P, S, Pc, dT, dTnext, Inj, Prod, FIM, Timers, Converged, CoarseGrid, Grid] = ...
                    FIMNonLinearSolver...
                (P0, S0, K, Trx, Try, Grid, Fluid, Inj, Prod, FIM, dT, Ndt, CoarseGrid, ADMSettings);
                FIM.MaxIter = maxiteration;
            else
                dT = min(dTnext, maxdT(index));
                % Newton-loop
                [P, S, Pc, dT, dTnext, Inj, Prod, FIM, Timers, Converged, CoarseGrid, Grid] = ...
                    FIMNonLinearSolver...
                (P0, S0, K, Trx, Try, Grid, Fluid, Inj, Prod, FIM, dT, Ndt, CoarseGrid, ADMSettings);
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
    Ndt = Ndt+1;
    
    %% COMPUTE OIL and WATER productions
    NwProduction(1, Ndt) =  t/(3600*24); 
    WProduction(1, Ndt) = t/(3600*24);
    for w=1:length(Prod)
        NwProduction(w+1, Ndt) = NwProduction(w+1, Ndt-1) - Prod(w).qnw*dT/(3600*24);
        WProduction(w+1,Ndt) = WProduction(w+1, Ndt-1) - Prod(w).qw*dT/(3600*24);
    end
    
    %% %%%%%%%%% POST PROCESSING and OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Print solution to a file at fixed intervals
    if (t == Tstops(index))
        disp(['Printing solution to file at  ' num2str((t)/(3600*24),4) ' days'])
        Saturations(:,index) = reshape(S, Grid.N, 1);
        Pressures(:,index) = reshape(P, Grid.N, 1);
    end   
    
    %%%%%%%%%%%%%%PLOT SOLUTION%%%%%%%%%%%%%
    switch (Options.PlotSolution)
        case('Matlab')
            if (t == Tstops(index))
                if ADMSettings.active
                    Plotting_ADM
                else
                    Plotting(Grid, P, Pc, S, Fluid, 'red', 'blue', Prod, Inj);
                end
                index = index +1;
            end
        case('VTK')
            if (t == Tstops(index))
                Write2VTK(Directory, Problem, vtkcount, Grid, K, P, S, Pc, ADMSettings.active, CoarseGrid, ADMSettings.maxLevel, 0);
                vtkcount = vtkcount + 1;
                index = index +1;
            end
    end
    
    %%%%%%%%%%%%%TImers of each timestep%%%%%%%%
    TimerTimestep(Ndt - 1) = toc(tstart);
    if (strcmp(Strategy, 'Sequential') == 1)
        TimerPressure(Ndt-1) = sum(Timers.ptimer);
        TimerBalance(Ndt-1) = sum(Timers.btimer);
        TimerSaturation(Ndt-1) = sum(Timers.stimer);
    else
        TimerConstruct(Ndt-1) = sum(Timers.Construct);
        TimerSolve(Ndt-1) = sum(Timers.Solve);
        if (ADMSettings.active == 1)
            TimerRP(Ndt-1) = Timers.RP;
        end
    end
end