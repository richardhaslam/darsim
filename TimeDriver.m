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
    TimerInner = zeros(TimeSteps,1);
    TimerRP = zeros(TimeSteps,1);
end

%%%%% START THE TIME LOOP %%%%%

t = 0;
Ndt = 1; 
Converged = 0;
index = 1;
Saturations = zeros(Grid.N, 10);
Pressures = zeros(Grid.N, 10);

Injection.time = zeros(1, TimeSteps);
Injection.Phase.W = zeros(length(Prod), TimeSteps);
Injection.Phase.Nw = zeros(length(Prod), TimeSteps);
Injection.Component.z1 = zeros(length(Prod), TimeSteps);
Injection.Component.z2 = zeros(length(Prod), TimeSteps);

Production.time = zeros(1, TimeSteps);
Production.Phase.W = zeros(length(Prod), TimeSteps);
Production.Phase.Nw = zeros(length(Prod), TimeSteps);
Production.Component.z1 = zeros(length(Prod), TimeSteps);
Production.Component.z2 = zeros(length(Prod), TimeSteps);
vtkcount = 1;
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
    switch (Strategy)
        case ('Sequential')
            disp('--------------Sequential non-linear solver-----');
            if ~Sequential.ImpSat
                Sequential.MaxExtIter = 1;
            end
            [Status, Pc, Inj, Prod, dT, Converged, Timers, Sequential.ImplicitSolver] =...
                SequentialStrategy(Status, K, Grid, Fluid, Inj, Prod, Sequential, Ndt, maxdT(index));
        case ('FIM')
            disp('------------FIM Non-linear solver--------------');
            FIM.timestep (Ndt) = Ndt;
            if (Ndt==1)
             
                %Plot initial conditions
                switch (Options.PlotSolution)
                    case('Matlab')
                        if ADMSettings.active
                            Plotting_ADM
                        else
                            Plotting(Grid, Status, Fluid, 'green', 'green', Inj, Prod);
                        end
                    case('VTK')
                        Write2VTK(Directory, Problem, vtkcount, Grid, K, Status, ADMSettings.active, CoarseGrid, ADMSettings.maxLevel, 1);
                        Wells2VTK(Grid, Inj, Prod, Directory, Problem);
                        vtkcount = vtkcount + 1;
                end
               
                % Use IMPES to estimate timestep size for the 1st timestep
                [~, U, ~, Wells, ~, ~] = PressureSolver(Grid, Inj, Prod, Fluid, Status.s, K);

                Grid.CFL = 1e-1; %Keep first timestep to be small
                maxiteration = FIM.MaxIter;
                FIM.MaxIter = 50;
                dTnext = timestepping(Fluid, Grid, U);
                dT = min(dTnext, maxdT(index));
                
                [Status, dT, dTnext, Inj, Prod, FIM, Timers, Converged, CoarseGrid, Grid] = ...
                    FIMNonLinearSolver...
                (Status, K, Grid, Fluid, Inj, Prod, FIM, dT, Ndt, CoarseGrid, ADMSettings, Directory, Problem, FlashSettings);
                FIM.MaxIter = maxiteration;
            else
                dT = min(dTnext, maxdT(index));

                %Non-linear solver
                [Status, dT, dTnext, Inj, Prod, FIM, Timers, Converged, CoarseGrid, Grid] = ...
                    FIMNonLinearSolver...
                (Status, K, Grid, Fluid, Inj, Prod, FIM, dT, Ndt, CoarseGrid, ADMSettings, Directory, Problem, FlashSettings);
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
    
    %% COMPUTE Phase and Components Injections and Production
    Injection.time(Ndt) = t/(3600*24);
    for w=1:length(Inj)
        Injection.Phase.W(w, Ndt) = Injection.Phase.W(w, Ndt-1) + Inj(w).qw*dT/(3600*24);
        Injection.Phase.Nw(w, Ndt) = Injection.Phase.Nw(w, Ndt-1) + Inj(w).qnw*dT/(3600*24);
        Injection.Component.z1(w, Ndt) = Injection.Component.z1(w, Ndt-1) + Inj(w).qz1*dT/(3600*24);
        Injection.Component.z2(w, Ndt) = Injection.Component.z2(w, Ndt-1) + Inj(w).qz2*dT/(3600*24);
    end
    
    Production.time(Ndt) = t/(3600*24);
    for w=1:length(Prod)
        Production.Phase.W(w, Ndt) = Production.Phase.W(w, Ndt-1) - Prod(w).qw*dT/(3600*24);
        Production.Phase.Nw(w, Ndt) = Production.Phase.Nw(w, Ndt-1) - Prod(w).qnw*dT/(3600*24);
        Production.Component.z1(w, Ndt) = Production.Component.z1(w, Ndt-1) - Prod(w).qz1*dT/(3600*24);
        Production.Component.z2(w, Ndt) = Production.Component.z2(w, Ndt-1) - Prod(w).qz2*dT/(3600*24);
    end
    
    
    %%%%%%%%%%%%%%PLOT SOLUTION%%%%%%%%%%%%%
    if (t == Tstops(index))
        disp(['Printing solution to file at  ' num2str((t)/(3600*24),4) ' days'])
        Saturations(:,index) = Status.s;
        Pressures(:,index) = Status.p;
        switch (Options.PlotSolution)
            case('Matlab')
                if ADMSettings.active
                    Plotting_ADM
                else
                    Plotting(Grid, Status, Fluid, 'red', 'green', Inj, Prod);
                end
            case('VTK')
                if (t == Tstops(index))
                    Write2VTK(Directory, Problem, vtkcount, Grid, K, Status, ADMSettings.active, CoarseGrid, ADMSettings.maxLevel, 0);
                    vtkcount = vtkcount + 1;
                end
        end
        index = index +1;
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
        TimerInner(Ndt-1) = sum(Timers.Inner);
        if (ADMSettings.active == 1)
            TimerRP(Ndt-1) = Timers.RP;
        end
    end
end