%FIM non-linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last Modified: 18 May 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIM non-linear solver
% solves the non-linear equations at timestep n with a fully implicit
% discretization.

%Input Variables:
%   P0: pressure at previous timestep
%   S0: saturation at previous time-step
%   K: permeability
%   Trx: rock transmissibility in x direction
%   Try: rock transmissibility in y direction
%   Grid: fine-scale grid information
%   Fluid: Fluid information
%   Inj: injection wells info
%   Prod: production wells info
%   FIM: FIM non-linear solver settings
%   dt: timestep size
%   Ndt: timestep number
%   CoarseGrid: ADM Coarse Grids
%   ADMSettings: settings for ADM

%Output variables
%   P: pressure at current timestep
%   S: saturation at current timestep
%   Pc: capillary pressure at current timestep
%   dt: timestep size (it can be modified for convergence issues)
%   dtnext: timestep size for next timestep
%   FIM: stats of FIM solver
%   Timers: timers of contruction and solve
%   Converged: 1 if non-linear convergence is achieved 0 otherwise
%   CoarseGrid: for ADM knows which coarse cells are active 
%   Grid: for ADM knows which fine cells are active

function [P, S, Pc, dt, dtnext, Inj, Prod, FIM, Timers, Converged, CoarseGrid, Grid] = ...
                    FIMNonLinearSolver...
                (P0, S0, K, Trx, Try, Grid, Fluid, Inj, Prod, FIM, dt, Ndt, CoarseGrid, ADMSettings)
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Grid.N;
pv = Grid.por*Grid.Volume;
ADM.active = 0;
Tol = FIM.Tol;
Kvector = reshape(K(1,:,:), N, 1);

% Initialise objects
Converged=0;
p0 = reshape(P0, N, 1);
s0 = reshape(S0, N, 1);
chops=0;
while (Converged==0 && chops<=10)
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
    [Pc, dPc] = ComputePc(s, Fluid, Kvector, Grid.por);
    %Define updwind operators
    [UpWindNw, Unw] = UpwindOperator(Grid, P, Trx, Try);
    [UpWindW, Uw] = UpwindOperator(Grid, P-reshape(Pc, Nx, Ny), Trx, Try);
    
    % Compute residual
    [Residual, ~, TMatrixNw, TMatrixW] = FIMResidual(p0, s0, p, s, Pc, pv, dt, Trx, Try, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, Kvector, N, Nx, Ny);

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
        J = BuildJacobian(Grid, Kvector, TMatrixNw, TMatrixW, p, Mw, Mnw, dMw, dMnw, Unw, Uw, dPc, dt, Inj, Prod, UpWindNw, UpWindW);
        TimerConstruct(itCount) = toc(start1);
       
        % 2. Solve full system at nu+1: J(nu)*Delta(nu+1) = -Residual(nu)
        start2 = tic;
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
        [Pc, dPc] = ComputePc(s, Fluid, Kvector, Grid.por); 
        % Define updwind
        [UpWindNw, Unw] = UpwindOperator(Grid, P, Trx, Try);
        [UpWindW, Uw] = UpwindOperator(Grid, P-reshape(Pc, Nx, Ny), Trx, Try);
        
        % 4. Compute residual 
        [Residual, qtot] = FIMResidual(p0, s0, p, s, Pc, pv, dt, Trx, Try, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, reshape(K(1,:,:), N, 1), N, Nx, Ny);

        % 5. Check convergence criteria
        Converged = NewtonConvergence(itCount, Residual, qtot, Delta(N+1:end), Tol, N, ADM);
        
        itCount = itCount+1;
    end
    if (Converged == 0)
        dt = dt/10;
        chops = chops + 1;
    end
end

%Choose next time-step size
if itCount < 5 
    dtnext = 2*dt;
elseif itCount > 8
    dtnext = dt/2;
else
    dtnext = dt;
end

% Reshape S before quitting
S = reshape(s,Nx,Ny);
Pc = reshape(Pc,Nx,Ny);

%Average saturation in Coarse Blocks
if ADM.active
    for x = 1:ADM.level
        S = Average(S, CoarseGrid(x), Grid);
    end
end

%Compute Nwetting and wetting phase fluxes for production curves
for i=1:length(Prod)
    c = Prod(i).cells;
    switch (Prod(i).type)
        case('RateConstrained')
            Prod(i).qw = sum (Mw(c)./(Mw(c)+Mnw(c)).*Prod(i).q);
            Prod(i).qnw = sum (Prod(i).q - Prod(i).qw(c));
        case('PressureConstrained')
            Prod(i).qnw =  sum(Mnw(c).* Prod(i).PI .* Kvector(c).* (Prod(i).p - p(c)))/(pv*Grid.N)*3600*24;
            Prod(i).qw =   sum(Mw(c).* Prod(i).PI .* Kvector(c) .* (Prod(i).p - p(c)))/(pv*Grid.N)*3600*24;
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