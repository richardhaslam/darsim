%FIM non-linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last modified: 21 March 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, S, Uo, dt, FIM, Timers, Converged, Inj, Prod, CoarseGrid, Grid] = ...
                    FIMNonLinearSolver...
                (P0, S0, K, Trx, Try, Grid, Fluid, Inj, Prod, FIM, dt, Ndt, CoarseGrid, ADMSettings)
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Grid.N;
pv = Grid.por*Grid.Volume;
ADM.active = 0;
Tol = FIM.Tol;

% Initialise objects
Converged=0;
p0 = reshape(P0, N, 1);
s0 = reshape(S0, N, 1);
chops=0;
while (Converged==0 && chops<=20)
    s=s0; % I start from solution at previous timestep.
    P = P0;
    p = p0;
    [Mw, Mo, dMw, dMo] = Mobilities(s, Fluid);
    %Oil
    [A, qo, Uo] = DivergenceMatrix(Grid, P, K, Trx, Try, Inj, Prod);
    Ro = -max(qo,0)*Inj(1).Mo -A*Mo + pv/dt*((1-s)-(1-s0));
    UpWindO = UpwindOperator(Grid, Uo);
    %Water
    %Pc = ComputePc(S); %Compute capillary pressure for all cells
    Pw = P; %- Pc;
    [A, qw, Uw] = DivergenceMatrix(Grid, Pw, K, Trx, Try, Inj, Prod);
    Rw = -max(qw,0)*Inj(1).Mw - A*Mw + pv/dt*(s-s0); %I am injecting water only
    Residual = [Ro; Rw];
    UpWindW = UpwindOperator(Grid, Uw);

    if (ADMSettings.active == 1 && chops == 0)
        tic
        % Choose where to coarsen and build DLGR grid
        [DLGRGrid, CoarseGrid, Grid] = AdaptGrid(Grid, CoarseGrid, S0, Ro, Rw, ADMSettings.maxLevel, ADMSettings.tol);
        % Construct R & P based on DLGR grid
        [ADM.Rest, ADM.Prolp, ADM.Prols] = ConstructOperators(Grid, CoarseGrid, DLGRGrid);
        ADM.level = DLGRGrid.level(end);
        ADM.active = 1;
        Timers.RP = toc;        
    end
    
    % NEWTON LOOP
    % Initialise Timers
    TimerConstruct = zeros(FIM.MaxIter,1);
    TimerSolve = zeros(FIM.MaxIter, 1);
    itCount = 1;
    while ((Converged==0) && (itCount < FIM.MaxIter))
        start1 = tic;
        % 1. Build Jacobian Matrix for nu+1: everything is computed at nu
        J = BuildJacobian(Grid, K, Trx, Try, P, Mw, Mo, dMw, dMo, ...
            Uo, Uw, dt, Inj, Prod, UpWindO, UpWindW);
        TimerConstruct(itCount) = toc(start1);
        
       
        % 2 Solve full system at nu+1: J(nu)*Delta(nu+1) = -Residual(nu)
        start2 = tic;
        Delta = LinearSolver(J, Residual, N, ADM);
        TimerSolve(itCount) = toc(start2);
      
        % 2.c Update solution
        p = p + Delta(1:N);
        s = s + Delta(N+1:2*N);
        s = min(s,1);
        s = max(s,0);
        P = reshape(p, Nx, Ny);
        
        % 3. Compute residuals at nu
        [Mw, Mo, dMw, dMo]=Mobilities(s, Fluid);
        %Oil
        [A, qo, Uo] = DivergenceMatrix(Grid, P, K, Trx, Try, Inj, Prod);
        Ro = -max(qo,0)*Inj(1).Mo -A*Mo + pv/dt*((1-s)-(1-s0));
        %Water
        %Pc = ComputePc(S); %Compute capillary pressure for all cells
        Pw = P; %- Pc;
        [A, qw, Uw] = DivergenceMatrix(Grid, Pw, K, Trx, Try, Inj, Prod);
        Rw = -max(qw,0)*Inj(1).Mw - A*Mw + pv/dt*(s-s0); %I am injecting water only
        Residual = [Ro; Rw];

        % 4. Check convergence criteria
        Converged = NewtonConvergence(Residual, Delta, Tol, N, ADM);
        
        itCount = itCount+1;
    end
    if (Converged == 0)
        dt = dt/2;
        chops = chops + 1;
    end
end

% Reshape S before quitting
S = reshape(s,Nx,Ny);

%Average saturation in Coarse Blocks
if (ADM.active == 1)
    for x = 1:ADM.level
        S = Average(S, CoarseGrid(x), Grid);
    end
end

%% Stats and timers 
FIM.Iter(Ndt) = itCount-1;
FIM.Chops(Ndt) = chops;
if ADMSettings.active == 1
    FIM.ActiveCells(Ndt, :) = DLGRGrid.N';
end
Timers.Construct = TimerConstruct;
Timers.Solve = TimerSolve;
end