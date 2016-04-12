%FIM non-linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last Modified: 6 April 2016
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
while (Converged==0 && chops<=0)
    if (chops > 0)
        disp('Maximum number of iterations was reached: time-step was chopped');
        disp('Restart Newton loop');
        disp('         ||Residual||  Sat. delta');
    end
    
    s = s0; % I start fRnwm solution at previous timestep.
    P = P0;
    p = p0;
    
    %Update fluid prowperties 
    [Mw, Mnw, dMw, dMnw] = Mobilities(s, Fluid);
    [Pc, dPc] = ComputePc(s, Fluid);
    %Define updwind operators
    [UpWindNw, Unw] = UpwindOperator(Grid, P, Trx, Try);
    [UpWindW, Uw] = UpwindOperator(Grid, P-reshape(Pc, Nx, Ny), Trx, Try);
    
    % Compute residual
    [Residual, TMatrixNw, TMatrixW] = FIMResidual(p0, s0, p, s, Pc, pv, dt, Trx, Try, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, reshape(K(1,:,:), N, 1), N, Nx, Ny);

    % Build ADM Grid and objects
    if (ADMSettings.active == 1 && chops == 0)
        tic
        % Choose where to coarsen and build ADM grid
        [ADMGrid, CoarseGrid, Grid] = AdaptGrid(Grid, CoarseGrid, S0, Residual(1:N), Residual(N+1:end), ADMSettings.maxLevel, ADMSettings.tol);
        % Construct R & P based on ADM grid
        [ADM.Rest, ADM.Prolp, ADM.Prols] = ConstructOperators(Grid, CoarseGrid, ADMGrid);
        ADM.level = ADMGrid.level(end);
        ADM.active = 1;
        Timers.RP = toc;        
    end
    
    % NEWTON LOOP
    % Initialise Timers
    TimerConstruct = zeros(FIM.MaxIter,1);
    TimerSolve = zeros(FIM.MaxIter, 1);
    itCount = 1;
    while ((Converged==0) && (itCount <= FIM.MaxIter))
              
        % 1. Build Jacobian Matrix for nu+1: everything is computed at nu
        start1 = tic;
        J = BuildJacobian(Grid, K, TMatrixNw, TMatrixW, p, Mw, Mnw, dMw, dMnw, Unw, Uw, Pc, dPc, dt, Inj, Prod, UpWindNw, UpWindW);
        TimerConstruct(itCount) = toc(start1);
       
        % 2. Solve full system at nu+1: J(nu)*Delta(nu+1) = -Residual(nu)
        start2 = tic;
        %[J, Residual] = ForceFixedValuesAtBoundary(J, Residual, N);
        Delta = LinearSolver(J, Residual, N, ADM);
        TimerSolve(itCount) = toc(start2);
      
        % 2.c Update solution
        p = p + Delta(1:N);
        s = s + Delta(N+1:2*N);
        s = min(s,1);
        s = max(s,0);
        P = reshape(p, Nx, Ny);
        
        % 3. Update fluid properties
        [Mw, Mnw, dMw, dMnw] = Mobilities(s, Fluid);
        [Pc, dPc] = ComputePc(s, Fluid); 
        % Define updwind
        [UpWindNw, Unw] = UpwindOperator(Grid, P, Trx, Try);
        [UpWindW, Uw] = UpwindOperator(Grid, P-reshape(Pc, Nx, Ny), Trx, Try);
        
        % 4. Compute residual 
        Residual = FIMResidual(p0, s0, p, s, Pc, pv, dt, Trx, Try, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, reshape(K(1,:,:), N, 1), N, Nx, Ny);
        %[J, Residual] = ForceFixedValuesAtBoundary(J, Residual, N);

        % 5. Check convergence criteria
        Converged = NewtonConvergence(itCount, Residual, Delta(N+1:end), Tol, N, ADM);
        
        itCount = itCount+1;
    end
    if (Converged == 0)
        dt = dt/10;
        chops = chops + 1;
    end
end

% Reshape S before quitting
S = reshape(s,Nx,Ny);

%Average saturation in Coarse Blocks
if ADM.active
    for x = 1:ADM.level
        S = Average(S, CoarseGrid(x), Grid);
    end
end

%% Stats and timers 
FIM.Iter(Ndt) = itCount-1;
FIM.Chops(Ndt) = chops;
if ADMSettings.active
    FIM.ActiveCells(Ndt, :) = ADMGrid.N';
end
Timers.Construct = TimerConstruct;
Timers.Solve = TimerSolve;
end