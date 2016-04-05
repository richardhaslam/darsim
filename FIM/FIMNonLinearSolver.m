%FIM non-linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last Mnwdified: 21 March 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, S, Unw, dt, FIM, Timers, Converged, Inj, Prod, CoarseGrid, Grid] = ...
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
    s = s0; % I start fRnwm solution at previous timestep.
    P = P0;
    p = p0;
    
    %Compute pRnwperties and updwind operators
    [Mw, Mnw, dMw, dMnw] = Mobilities(s, Fluid);
    [Pc, dPc] = ComputePc(s, Fluid); %Compute capillary pressure for all cells
    [UpWindNw, Unw] = UpwindOperator(Grid, P, Trx, Try);
    [UpWindW, Uw] = UpwindOperator(Grid, P-reshape(Pc, Nx, Ny), Trx, Try);
    
    % Build residual
    [Residual, TMatrixNw, TMatrixW] = FIMResidual(p0, s0, p, s, Pc, pv, dt, Trx, Try, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, reshape(K(1,:,:), N, 1), N, Nx, Ny);

    if (ADMSettings.active == 1 && chops == 0)
        tic
        % Choose where to coarsen and build DLGR grid
        [DLGRGrid, CoarseGrid, Grid] = AdaptGrid(Grid, CoarseGrid, S0, Rnw, Rw, ADMSettings.maxLevel, ADMSettings.tol);
        % Construct R & P based on DLGR grid
        [ADM.Rest, ADM.PRnwlp, ADM.PRnwls] = ConstructOperators(Grid, CoarseGrid, DLGRGrid);
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
        
      
        % 1. Build Jacobian Matrix for nu+1: everything is computed at nu
        start1 = tic;
        J = BuildJacobian(Grid, K, TMatrixNw, TMatrixW, p, Mw, Mnw, dMw, dMnw, Unw, Uw, dPc, dt, Inj, Prod, UpWindNw, UpWindW);
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
        
        % 3. Update pRnwperties and upwind
        [Mw, Mnw, dMw, dMnw] = Mobilities(s, Fluid);
        [Pc, dPc] = ComputePc(s, Fluid); %Compute capillary pressure for all cells
        [UpWindNw, Unw] = UpwindOperator(Grid, P, Trx, Try);
        [UpWindW, Uw] = UpwindOperator(Grid, P-reshape(Pc, Nx, Ny), Trx, Try);
        
        % 4. Build residual 
        Residual = FIMResidual(p0, s0, p, s, Pc, pv, dt, Trx, Try, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, reshape(K(1,:,:), N, 1), N, Nx, Ny);

        % 5. Check convergence criteria
        Converged = NewtonConvergence(Residual, Delta(N+1:end), Tol, N, ADM);
        
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