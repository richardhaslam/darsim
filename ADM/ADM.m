%           ADM CODE                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Timers
TimerTimestep=zeros(TimeSteps,1);
TimerConstruct=zeros(TimeSteps,1);
TimerRP=zeros(TimeSteps,1);
TimerSolve=zeros(TimeSteps,1);

%% Construct Coarse Grids
CoarseGridding;

%% %%%START THE TIME LOOP%%%%%
Tstops = linspace(T/10, T, 10);
index = 1;
Saturations = zeros(FineGrid.N, 10);
Pressures = zeros(FineGrid.N, 10);
while (t<T && Ndt <= TimeSteps)
    tstart = tic;
    S0 = S;
    P0 = P;
    %Reset all fine cells to active
    FineGrid.Active = ones(FineGrid.N, 1);
    
    %Solve Timestep
    maxdT=Tstops - t;
    switch (Strategy)
        case ('Sequential')
        disp('No DLGR has been implemented for Sequential Simulations yet');
        break;
        case ('FIM')
            if (Ndt==1)
                disp('******************');
                disp('FIM strategy');
                disp(char(10));
            end
            FIM.timestep (Ndt) = Ndt;
            if (Ndt==1) % I only need to do this at the 1st timestep
                % Use IMPES as intial guess for pressure for the 1st timestep
                % Effective permeability
                [Mw, Mo]=Mobilities(S, Fluid);
                Mt=Mw+Mo;   %total mobility
                Kt=zeros(2, FineGrid.Nx, FineGrid.Ny);
                Kt(1,:,:)=reshape(Mt, 1, FineGrid.Nx, FineGrid.Ny).*K(1,:,:);   % x-direction
                Kt(2,:,:)=reshape(Mt, 1, FineGrid.Nx, FineGrid.Ny).*K(2,:,:);   % y-direction
                % PressureSolver as initial guess
                [P0, U, Wells, Ap, Ab, q] = PressureSolver(FineGrid, Kt, Inj, Prod);
                %[P0, U_mmsfv] = MMsFVPressureSolver(FineGrid, Inj, Prod, Kt, CoarseGrid, maxLevel);
                
                
                %Keep first timestep to be small
                FineGrid.CFL = 0.25/8;
                dT = timestepping(Fluid, S, FineGrid, U, Wells);
                %Compute timestep size based on CFL
                FineGrid.CFL = FIM.CFL;
                dT_CFL = timestepping(Fluid, S, FineGrid, U, Wells);
                
                %Compute rock transmissibility
                [Trx, Try] = ComputeTransmissibility(FineGrid, K);
                
                %Pressure Interpolator
                disp('Pressure interpolator - start computation');
                CoarseGrid = PressureInterpolator(FineGrid, Kt, CoarseGrid, maxLevel, ADMSettings.Pressure_Interpolator);
                disp('Pressure interpolator - end')
                disp(char(2));
                               
                %Solve first timestep
                %P0 = zeros(FineGrid.Nx, FineGrid.Ny);
                [P, S, U, dT, FIM, Timers, ActiveFine, CoarseGrid, FineGrid, Converged, Inj, Prod] = ...
                FullyImplicit_DLGR(P0, S0, K, Trx, Try, FineGrid, ...
                CoarseGrid, Fluid, Inj, Prod, FIM, dT, Options, Ndt, ADMSettings);
            else
                dT = min(dT_CFL, maxdT(index));
                [P, S, U, dT, FIM, Timers, ActiveFine, CoarseGrid, FineGrid, Converged, Inj, Prod, ADMProl] = ...
                    FullyImplicit_DLGR(P0, S0, K, Trx, Try, FineGrid, ...
                    CoarseGrid, Fluid, Inj, Prod, FIM, dT, Options, Ndt, ADMSettings);
            end
            
            %Average Solution in Coarse Blocks
            for x = 1:maxLevel
                [P_, S] = Average(P, S, CoarseGrid(x), FineGrid);
            end
    end
    
    %Check for convergence at the end of the timestep
    if (Converged==0) 
        switch(Strategy) 
            case('Sequential')
                disp(['The solution has not converged at timestep ' num2str(Ndt)]);
                break
            case('FIM')
                disp(['The solution has not converged at timestep '...
                    num2str(Ndt) ' even with 20 timestep chops']);
                break
        end
    end
    
    TimerTimestep(Ndt) = toc(tstart);
    
    %%%%%Increase time and timestep counter
    disp(['Timestep ' num2str(Ndt)]);
    disp(['Initial time: ' num2str(t/(3600*24),4) ' days -- Final time: ' num2str((t+dT)/(3600*24), 4) ' days, dT= ' num2str(dT) ' s']);
    disp(char(5));  
    t = t+dT;
    Ndt = Ndt+1;
    CumulativeTime(Ndt) = t/(3600*24);
    %Print solution to a file at fixed intervals
    if (t == Tstops(index))
        disp(['Printing solution to file at  ' num2str((t)/(3600*24),4) ' days'])
        Saturations(:,index) = reshape(S, FineGrid.N, 1);
        Pressures(:,index) = reshape(P, FineGrid.N, 1);
        index = index + 1;
    end
    
     %%%%%%%%%%%%%%PLOT SOLUTION%%%%%%%%%%%%%
    switch (Options.PlotSolution)
        case('Matlab')
            if (Grid.Nx == 1 || Grid.Ny == 1)
            Options.problem_1D=1;
            end
        Plotting_DLGR;
        case('VTK')
            Write2VTK(Directory, Problem, Ndt, Grid, K, P, S);
    end
    
    %%%%%%%%%%%%%Timers
    TimerRP(Ndt-1) = Timers.RP;
    TimerConstruct(Ndt-1) = sum(Timers.Construct);
    TimerSolve(Ndt-1) = sum(Timers.Solve);
end
%% %%%Output%%%%%%%%%%%%%%%%%
OutputStat_DLGR

%% Injection and Production data
%PlotInjProdCurves