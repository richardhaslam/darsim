%        Time loop driver             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last modified: 9 April 2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Timers and variables for statistics
TimerTimestep = zeros(TimeSteps,1);
CumulativeTime = zeros(TimeSteps, 1);
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
Tstops = linspace(T/10, T, 10);
index = 1;
Saturations = zeros(Grid.N, 10);
Pressures = zeros(Grid.N, 10);
vtkcount = 1;

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
            [P, S, Pc, dT, Converged, Timers, Sequential.ImplicitSolver] =...
                SequentialStrategy(S0, K, Grid, Fluid, Inj, Prod, Sequential, Ndt, maxdT(index));
        case ('FIM')
            disp('------------FIM Non-linear solver--------------');
            disp('        ||Residual||   Sat. delta');
            FIM.timestep (Ndt) = Ndt;
            if (Ndt==1)
                % Use IMPES as intial guess for pressure for the 1st timestep
                [~, U, Pc, Wells] = PressureSolver(Grid, Inj, Prod, Fluid, S, K);
                %[Pms, ~] = MMsFVPressureSolver(Grid, Inj, Prod, K, Fluid, S, CoarseGrid, maxLevel);
                
                %Keep first timestep to be small
                Grid.CFL = 0.25/8;
                dT = timestepping(Fluid, S, Grid, U, Wells);
               
                %Compute timestep size based on CFL
                Grid.CFL = FIM.CFL;
                dT_CFL = timestepping(Fluid, S, Grid, U, Wells);
               
                %Compute rock transmissibility
                [Trx, Try] = ComputeTransmissibility(Grid, K);
                
                %Non-linear solver
                [P, S, Pc, U, dT, FIM, Timers, Converged, Inj, Prod, CoarseGrid, Grid] = ...
                    FIMNonLinearSolver...
                (P0, S0, K, Trx, Try, Grid, Fluid, Inj, Prod, FIM, dT, Ndt, CoarseGrid, ADMSettings);
            else
                dT = min(dT_CFL, maxdT(index));
                % Newton-loop
                [P, S, Pc, U, dT, FIM, Timers, Converged, Inj, Prod, CoarseGrid, Grid] = ...
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
    CumulativeTime(Ndt) = t/(3600*24);
    
    %% %%%%%%%%% POST PROCESSING and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Print solution to a file at fixed intervals
    if (t == Tstops(index))
        disp(['Printing solution to file at  ' num2str((t)/(3600*24),4) ' days'])
        Saturations(:,index) = reshape(S, Grid.N, 1);
        Pressures(:,index) = reshape(P, Grid.N, 1);
        index = index + 1;
    end
      
    %%%%%%%%%%%%%%PLOT SOLUTION%%%%%%%%%%%%%
    switch (Options.PlotSolution)
        case('Matlab')
            if ADMSettings.active
                Plotting_ADM
            else
                Plotting;
            end
        case('VTK')
            if (mod(Ndt,5)==0)
                Write2VTK(Directory, Problem, vtkcount, Grid, K, P, S, Pc, ADMSettings.active, CoarseGrid, maxLevel);
                vtkcount = vtkcount + 1;
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