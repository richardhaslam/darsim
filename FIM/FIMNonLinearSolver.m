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
                (Status0, K, Grid, Fluid, Inj, Prod, FIM, dt, Ndt, CoarseGrid, ADMSettings, Directory, Problem)
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Grid.N;
pv = Grid.por*Grid.Volume;
ADM.active = 0;
Tol = FIM.Tol;
Kvector = reshape(K(1,:,:), N, 1);

% Initialise objects
Converged=0;
chops=0;
while (Converged==0 && chops<=10)
    if (chops > 0)
        disp('Maximum number of iterations was reached: time-step was chopped');
        disp('Restart Newton loop');
        disp('         ||Residual||  Sat. delta');
    end
    
    Status.p = Status0.p;
    Status.s = Status0.s; % I start from solution at previous timestep.
    Status.x1 = Status0.x1;
    Status.z = Status0.z;
    P = reshape(Status.p,Grid.Nx,Grid.Ny);
     
    %Update fluid prowperties 
    [Mw, Mnw, dMw, dMnw] = Mobilities(Status.s, Fluid);
    [Rho, dRho] = LinearDensity(Status.p, Fluid.c, Fluid.rho); 
    [Status.pc, dPc] = ComputePc(Status.s, Fluid, Kvector, Grid.por);
    %Define updwind operators
    [UpWindW, Uw] = UpwindOperator(Grid, P-reshape(Status.pc, Nx, Ny));
    [UpWindNw, Unw] = UpwindOperator(Grid, P);
    
    % Compute residual
    [Residual1, TMatrixNw, TMatrixW] = FIMResidual(Status0, Status, Grid, dt, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, Kvector);
    [Residual, TMatrix1, TMatrix2] = FIMResidualComp(Status0, Status, dt, Grid, Kvector, Fluid, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod);
    
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
        J1 = BuildJacobian(Grid, Kvector, TMatrixNw, TMatrixW, Status.p, Mw, Mnw, dMw, dMnw, Unw, Uw, dPc, dt, Inj, Prod, UpWindNw, UpWindW);
        J = BuildJacobianComp(Grid, Kvector, TMatrix1, TMatrix2, Status, Mw, Mnw, dMw, dMnw, Rho, dRho, Uw, Unw, dPc, dt, Inj, Prod, UpWindW, UpWindNw);
        %max(Residual - Residual1)
        %J - J1
        TimerConstruct(itCount) = toc(start1);
       
        % 2. Solve full system at nu+1: J(nu)*Delta(nu+1) = -Residual(nu)
        start2 = tic;
        [Delta, Delta_c] = LinearSolver(J, Residual, N, ADM);
        TimerSolve(itCount) = toc(start2);
      
        % 2.c Update solution
        Status.p = Status.p + Delta(1:N);
        Status.s = Status.s + Delta(N+1:2*N);
        Status.s = min(Status.s,1);
        Status.s = max(Status.s,0);
        P = reshape(Status.p, Nx, Ny);
        
        % 3. Update fluid properties
        [Mw, Mnw, dMw, dMnw] = Mobilities(Status.s, Fluid);
        [Status.pc, dPc] = ComputePc(Status.s, Fluid, Kvector, Grid.por); 
        % Define updwind
        [UpWindNw, Unw] = UpwindOperator(Grid, P);
        [UpWindW, Uw] = UpwindOperator(Grid, P-reshape(Status.pc, Nx, Ny));
        
        % Inner Update and Flash
        %[Status] = Inner_Update(Status,Fluid,FlashSettings,Grid);
          
        %Print residual if required
%         if (Ndt == 50000 && ADMSettings.active == 1)
%             Residualc = RestrictResidual(Residual, ADM.Rest, Grid.N, ADM.level);
%             Residualc = Prolong(Residualc, ADM.Prols, ADM.Prols, ADM.level);
%             WriteADMResiduals2VTK(Directory, Problem, itCount, Grid, abs(Residual), abs(Residualc), CoarseGrid, ADM.level)
%         end
        
        % 4. Compute residual 
        [Residual1, TMatrixNw, TMatrixW] = FIMResidual(Status0, Status, Grid, dt, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, Kvector);
        [Residual, TMatrix1, TMatrix2] = FIMResidualComp(Status0, Status, dt, Grid, Kvector, Fluid, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod);
        
        % 5. Check convergence criteria
        Converged = NewtonConvergence(itCount, Residual, Delta, Status.p, Tol, N, ADM, Delta_c);
        
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
            Prod(i).qnw =  sum(Mnw(c).* Prod(i).PI .* Kvector(c).* (Prod(i).p - Status.p(c)))/(pv*Grid.N)*3600*24;
            Prod(i).qw =   sum(Mw(c).* Prod(i).PI .* Kvector(c) .* (Prod(i).p - Status.p(c)))/(pv*Grid.N)*3600*24;
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