function [P, S, U, dt, FIM, Timers, ActiveFine, CoarseGrid, FineGrid, Converged, Inj, Prod, Prolp] = ...
    FullyImplicit_DLGR(P0, S0, K, Trx, Try, FineGrid, CoarseGrid,...
    Fluid, Inj, Prod, FIM, dt, Options, Ndt, ADMSettings)
%% DLGR: FULLY IMPLICIT STRATEGY
Nx = FineGrid.Nx;
Ny = FineGrid.Ny;
N = FineGrid.N;
pv = FineGrid.por*FineGrid.Volume;
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
    [A, q, U] = DivergenceMatrix(FineGrid, P, K, Trx, Try, Inj, Prod);
    Ro = -A*Mo + pv/dt*((1-s)-(1-s0));
    Rw = -max(q,0)*Inj.Mw - A*Mw + pv/dt*(s-s0); %I am injecting water only
    Residual = [Ro; Rw];
    UpWind = UpwindOperator(FineGrid, U);

    if (chops == 0)
        tic
        % Choose where to coarsen and build DLGR grid
        [DLGRGrid, ActiveFine, CoarseGrid, FineGrid] = AdaptGrid(FineGrid, CoarseGrid, S0, Ro, Rw, ADMSettings.maxLevel, ADMSettings.tol);
        % Construct R & P based on DLGR grid
        [Rest, Prolp, Prols] = ConstructOperators(FineGrid, CoarseGrid, DLGRGrid);       
        Timers.RP = toc;        
                
        % Plot Residuals at beginning of timestep
        if (Options.PlotResiduals == 1)
            PlotResiduals(Ro, Rw, FineGrid);
        end
    end
    
    % Timers
    TimerConstruct = zeros(FIM.MaxIter,1);
    TimerSolve = zeros(FIM.MaxIter, 1);
    
    % NEWTON LOOP
    itCount = 1;
    while ((Converged==0) && (itCount < FIM.MaxIter))
        tic
        % 1. Build Jacobian Matrix for nu+1: everything is computed at nu
        J = BuildJacobian(FineGrid, K, Trx, Try, P, Mw, Mo, dMw, dMo, ...
            U, dt, Inj, Prod, UpWind);
        TimerConstruct(itCount) = toc;
        
        tic
        % 2.a Solve full system at nu+1: J(nu)*Delta(nu+1) = -Residual(nu)
        [J_c, Residual_c] = Restrict(J, Residual, Rest, Prolp, Prols, N, DLGRGrid.level(end));
        Delta_c = -J_c\Residual_c;
        TimerSolve(itCount) = toc;
        
        % 2.b Prolong Solution
        Delta = Prolong(Delta_c, Prolp, Prols, DLGRGrid.level(end));
        
        % 2.c Update solution
        p = p + Delta(1:N);
        s = s + Delta(N+1:2*N);
        s = min(s,1);
        s = max(s,0);
        P = reshape(p, Nx, Ny);
        
        % 3. Compute residuals at nu
        [Mw, Mo, dMw, dMo]=Mobilities(s, Fluid);
        [A, q, U] = DivergenceMatrix(FineGrid, P, K, Trx, Try, Inj, Prod);
        Ro = -A*Mo + pv/dt*((1-s)-(1-s0));
        Rw = -max(q,0)*Inj.Mw - A*Mw + pv/dt*(s-s0); %I am injecting water only
        Residual = [Ro; Rw];
        
        % 4. Check convergence criteria
        Norm1 = norm(Residual_c, inf);
        %Norm2 = norm(Delta_c(DLGRGrid.N:2*DLGRGrid.N), inf);
        if (Norm1 < Tol) %&& Norm2 < Tol)
            Converged = 1;
        end
        itCount = itCount+1;
    end
    if (Converged == 0)
        dt = round(dt/2);
        chops = chops + 1;
    end
end

% Reshape S before quitting
S = reshape(s,Nx,Ny);
% Update production and injection data in wells
Inj.water(Ndt+1) = Inj.water(Ndt) + dt*Inj.PI*(K(1,Inj.x,Inj.y)...
    *K(2, Inj.x, Inj.y))^0.5*(Inj.p-P(Inj.x,Inj.y))*Inj.Mw;
Prod.water(Ndt+1) = Prod.water(Ndt) - dt*Prod.PI*(K(1,Prod.x,Prod.y)...
    *K(2, Prod.x, Prod.y))^0.5*(Prod.p-P(Prod.x,Prod.y))*Mw((Prod.y-1)*Nx+Prod.x);
Prod.oil(Ndt+1) =  Prod.oil(Ndt) - dt*Prod.PI*(K(1,Prod.x,Prod.y)...
    *K(2, Prod.x, Prod.y))^0.5*(Prod.p-P(Prod.x,Prod.y))*Mo((Prod.y-1)*Nx+Prod.x);
%% Stats
FIM.Iter(Ndt) = itCount-1;
FIM.Chops(Ndt) = chops;
FIM.ActiveCells(Ndt, :) = DLGRGrid.N';
Timers.Construct = TimerConstruct;
Timers.Solve = TimerSolve;
end