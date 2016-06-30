%FIM non-linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini and Barnaby Fryer
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
%   Status
%   dt: timestep size (it can be modified for convergence issues)
%   dtnext: timestep size for next timestep
%   FIM: stats of FIM solver
%   Timers: timers of contruction and solve
%   Converged: 1 if non-linear convergence is achieved 0 otherwise
%   CoarseGrid: for ADM knows which coarse cells are active 
%   Grid: for ADM knows which fine cells are active

function [Status, dt, dtnext, Inj, Prod, FIM, Timers, Converged, CoarseGrid, Grid] = ...
                    FIMNonLinearSolver...
                (Status0, K, Grid, Fluid, Inj, Prod, FIM, dt, Ndt, CoarseGrid, ADMSettings, Directory, Problem, FlashSettings)
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Grid.N;
ADM.active = 0;
Tol = FIM.Tol;
Kvector = reshape(K(1,:,:), N, 1);

% Initialise objects
Converged=0;
CompConverged = 0;
chops=0;
while ((Converged==0 || CompConverged == 0) && chops<=10)
    Status.p = Status0.p;
    Status.s = Status0.s; % I start from solution at previous timestep.
    Status.x1 = Status0.x1;
    Status.z = Status0.z;
    P = reshape(Status.p,Grid.Nx,Grid.Ny);
     
    %Update fluid prowperties 
    [Mw, Mnw, dMw, dMnw] = Mobilities(Status.s, Fluid);
    [Status.rho, dRho] = LinearDensity(Status.p, Fluid.c, Fluid.rho); 
    [Status.pc, dPc] = ComputePc(Status.s, Fluid, Kvector, Grid.por);
    %Define updwind operators
    [UpWindW, Uw] = UpwindOperator(Grid, P-reshape(Status.pc, Nx, Ny));
    [UpWindNw, Unw] = UpwindOperator(Grid, P);
    
    % Compute residual
    %[Residual1, TMatrixNw, TMatrixW] = FIMResidual(Status0, Status, Grid, dt, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, Kvector);
    [Residual, TMatrix1, TMatrix2, TMatrixW] = FIMResidualComp(Status0, Status, dt, Grid, Kvector, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod);
    
    %Print some info to the screen 
    if (chops > 0)
        disp('Maximum number of iterations was reached: time-step was chopped');
        disp(['Restart Newton loop dt = ', num2str(dt)]);
    end
    disp(['Initial residual norm: ', num2str(norm(Residual, inf))]);
    disp('');
    disp('        ||Residual||   ||delta p||   ||delta S||');
    
    % Build ADM Grid and objects
    if (ADMSettings.active == 1 && chops == 0)
        tic
        % Choose where to coarsen and build ADM grid
        [ADMGrid, CoarseGrid, Grid] = AdaptGrid(Grid, CoarseGrid, reshape(Status0.s,Nx,Ny), Residual(1:N), Residual(N+1:end), ADMSettings.maxLevel, ADMSettings.tol);
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
    TimerInner = zeros(FIM.MaxIter, 1);
    itCount = 1;
    while ((Converged==0 || CompConverged == 0)  && (itCount <= FIM.MaxIter))
              
        % 1. Build Jacobian Matrix for nu+1: everything is computed at nu
        start1 = tic;
        %J1 = BuildJacobian(Grid, Kvector, TMatrixNw, TMatrixW, Status.p, Mw, Mnw, dMw, dMnw, Unw, Uw, dPc, dt, Inj, Prod, UpWindNw, UpWindW);
        J = BuildJacobianComp(Grid, Kvector, TMatrix1, TMatrix2, TMatrixW, Status, Mw, Mnw, dMw, dMnw, dRho, Uw, Unw, dPc, dt, Inj, Prod, UpWindW, UpWindNw);
        TimerConstruct(itCount) = toc(start1);
       
        % 2. Solve full system at nu+1: J(nu)*Delta(nu+1) = -Residual(nu)
        start2 = tic;
        [Delta, Delta_c] = LinearSolver(J, Residual, N, ADM);
        TimerSolve(itCount) = toc(start2);
      
        % 2.c Update solution
        Status.p = Status.p + Delta(1:N);
        Status.s = Status.s + Delta(N+1:2*N);
        Status.p = max(Status.p,0);
        
        P = reshape(Status.p, Nx, Ny);
     
        % 2.d Update solution based on phase split
        start3 = tic;
        [Status.rho, dRho] = LinearDensity(Status.p, Fluid.c, Fluid.rho);
        [Status, CompConverged] = CompositionUpdate(Status, Fluid, Grid, FlashSettings);
        TimerInner(itCount) = toc(start3);
        
        % 3. Update fluid properties
        [Mw, Mnw, dMw, dMnw] = Mobilities(Status.s, Fluid);
        [Status.pc, dPc] = ComputePc(Status.s, Fluid, Kvector, Grid.por); 
        % Define updwind
        [UpWindNw, Unw] = UpwindOperator(Grid, P);
        [UpWindW, Uw] = UpwindOperator(Grid, P-reshape(Status.pc, Nx, Ny));
        
        % 4. Compute residual 
        %[Residual1, TMatrixNw, TMatrixW] = FIMResidual(Status0, Status, Grid, dt, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, Kvector);
        [Residual, TMatrix1, TMatrix2, TMatrixW] = FIMResidualComp(Status0, Status, dt, Grid, Kvector, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod);
        
        % 5. Check convergence criteria
        Converged = NewtonConvergence(itCount, Residual, Delta, Status.p, Tol, N, ADM, Delta_c);
        itCount = itCount+1;
    end
    if (Converged == 0 || CompConverged == 0)
        dt = dt/2;
        chops = chops + 1;
    end
end

%Choose next time-step size
if itCount < 12 
    dtnext = 2*dt;
elseif itCount > 15
    dtnext = dt;
else
    dtnext = dt;
end

%Average saturation in Coarse Blocks
if ADM.active
    for x = 1:ADM.level
        Status.s = Average(Status.s, CoarseGrid(x), Grid);
    end
    Status.s = reshape(Status.s,Nx*Ny,1);
end

%Compute Injection and production fluxes for Injection  and Production curves
for i=1:length(Inj)
    c = Inj(i).cells;
    switch (Inj(i).type)
        case('RateConstrained')
            Inj(i).qw = sum (Mw(c)./(Mw(c)+Mnw(c)).*Inj(i).q);
            Inj(i).qnw = sum (Inj(i).q - Inj(i).qw(c));
        case('PressureConstrained')
            %Phases
            Inj(i).qw =   sum(Mw(c).* Inj(i).rho(c,1).* Inj(i).PI .* Kvector(c) .* (Inj(i).p - Status.p(c)))*3600*24;
            Inj(i).qnw =  sum(Mnw(c).* Inj(i).rho(c,2) .* Inj(i).PI .* Kvector(c).* (Inj(i).p - Status.p(c)))*3600*24;
            %Components
            Inj(i).qz1 = sum(Inj(i).x1(1) .* Inj(i).Mw .* Inj.rho(1) .* Inj(i).PI .* Kvector(c) .* (Inj(i).p - Status.p(c)))*3600*24 +... 
                          sum(Inj(i).x1(2) .* Inj(i).Mo(c) .* Inj(i).rho(2) .* Inj(i).PI .* Kvector(c) .* (Inj(i).p - Status.p(c)))*3600*24;
            Inj(i).qz2 = sum((1 - Inj(i).x1(1)) .* Inj(i).Mw .* Inj(i).rho(1) .* Inj(i).PI .* Kvector(c) .* (Inj(i).p - Status.p(c)))*3600*24 +... 
                          sum((1 - Inj(i).x1(2)) .* Inj(i).Mo .* Inj(i).rho(2) .* Inj(i).PI .* Kvector(c) .* (Inj(i).p - Status.p(c)))*3600*24; 
    end
end

for i=1:length(Prod)
    c = Prod(i).cells;
    switch (Prod(i).type)
        case('RateConstrained')
            Prod(i).qw = sum (Mw(c)./(Mw(c)+Mnw(c)).*Prod(i).q);
            Prod(i).qnw = sum (Prod(i).q - Prod(i).qw(c));
        case('PressureConstrained')
            %Phases
            Prod(i).qw =   sum(Mw(c).* Status.rho(c,1).* Prod(i).PI .* Kvector(c) .* (Prod(i).p - Status.p(c)))*3600*24;
            Prod(i).qnw =  sum(Mnw(c).* Status.rho(c,2) .* Prod(i).PI .* Kvector(c).* (Prod(i).p - Status.p(c)))*3600*24;
            %Components
            Prod(i).qz1 = sum(Status.x1(c, 1) .* Mw(c) .* Status.rho(c,1) .* Prod(i).PI .* Kvector(c) .* (Prod(i).p - Status.p(c)))*3600*24 +... 
                          sum(Status.x1(c, 2) .* Mnw(c) .* Status.rho(c,2) .* Prod(i).PI .* Kvector(c) .* (Prod(i).p - Status.p(c)))*3600*24;
            Prod(i).qz2 = sum((1 - Status.x1(c, 1)) .* Mw(c) .* Status.rho(c,1) .* Prod(i).PI .* Kvector(c) .* (Prod(i).p - Status.p(c)))*3600*24 +... 
                          sum((1 - Status.x1(c, 2)) .* Mnw(c) .* Status.rho(c,2) .* Prod(i).PI .* Kvector(c) .* (Prod(i).p - Status.p(c)))*3600*24; 
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
Timers.Inner = TimerInner;
end