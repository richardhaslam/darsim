%FIM non-linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, S, Uo, dt, FIM, Timers, Converged, Inj, Prod] = ...
    FullyImplicit(P0, S0, K, Trx, Try, Grid, Fluid, Inj, Prod, FIM, dt, Options, Ndt)

%FULLY IMPLICIT STRATEGY
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Grid.N;
pv = Grid.por*Grid.Volume;
Tol = FIM.Tol;

%Initialise objects
Converged=0;
p0 = reshape(P0,N,1);
s0 = reshape(S0, N, 1);
chops=0;
while (Converged==0 && chops <= 20)
    s=s0; % I start from solution at previous timestep.
    P = P0;
    p = p0;
    [Mw, Mo, dMw, dMo] = Mobilities(s, Fluid);
    %Oil
    [A, qo, Uo] = DivergenceMatrix(Grid, P, K, Trx, Try, Inj, Prod);
    Ro = -max(qo,0)*Inj.Mo -A*Mo + pv/dt*((1-s)-(1-s0));
    UpWindO = UpwindOperator(Grid, Uo);
    %Water
    %Pc = ComputePc(S); %Compute capillary pressure for all cells
    Pw = P; %- Pc;
    [A, qw, Uw] = DivergenceMatrix(Grid, Pw, K, Trx, Try, Inj, Prod);
    Rw = -max(qw,0)*Inj.Mw - A*Mw + pv/dt*(s-s0); %I am injecting water only
    Residual = [Ro; Rw];
    UpWindW = UpwindOperator(Grid, Uw);
    
    %Plot Residuals at beginning of timestep
    if (Options.PlotResiduals==1) 
        PlotResiduals(Ro, Rw, Grid);
    end
    
    %Timers
    TimerConstruct = zeros(FIM.MaxIter,1);
    TimerSolve = zeros(FIM.MaxIter, 1);
    
    itCount = 1;
    while ((Converged==0) && (itCount < FIM.MaxIter))
        tic
        % 1. Build Jacobian Matrix for nu+1: everything is computed at nu
        J = BuildJacobian(Grid, K, Trx, Try, P, Mw, Mo, dMw, dMo,...
            Uo, Uw, dt, Inj, Prod, UpWindO, UpWindW);
        TimerConstruct(itCount) = toc;
        %Rescale for precision purposes
        J = 1e6 * J;
        Residual = 1e6 * Residual;
        %CheckEllipticSystem(J(1:Grid.N, 1:Grid.N), Grid.Nx);
        TimerConstruct(itCount) = toc;
       
        tic
        % 2.a Solve full system at nu+1: J(nu)*Delta(nu+1) = -Residual(nu)
        Delta = -J\Residual;
        %Delta = gmres(J, Residual, 100,1e-3,100); Does not seem to work
        TimerSolve(itCount) = toc;
        
        % 2.b Update solution
        p = p + Delta(1:N);
        s = s + Delta(N+1:2*N);
        s = min(s,1-Fluid.swc);
        s = max(s,Fluid.swc);        
        P = reshape(p,Nx,Ny);
        
        % 3. Compute residuals at nu
        [Mw, Mo, dMw, dMo]=Mobilities(s, Fluid);
        %Oil
        [Ao, qo, Uo] = DivergenceMatrix(Grid, P, K, Trx, Try, Inj, Prod);
        Ro = -max(qo,0)*Inj.Mo - Ao*Mo + pv/dt*((1-s)-(1-s0));
        %Water
        %Pc = ComputePc(S); %Compute capillary pressure for all cells
        Pw = P; %- Pc;
        [Aw, qw, Uw] = DivergenceMatrix(Grid, Pw, K, Trx, Try, Inj, Prod);
        Rw = -max(qw,0)*Inj.Mw - Aw*Mw + pv/dt*(s-s0); %I am injecting water only
        Residual = [Ro; Rw];
        
        % 4. Check convergence criteria
        Norm1 = max(norm(Ro, inf), norm(Rw, inf));
        Norm2 = norm(Delta(N+1:2*N), inf);
        Norm3 = norm(Delta(1:N)/Inj.p,inf);
        if (Norm1 < Tol && Norm2 < Tol && Norm3 < Tol)
            Converged = 1;
        end
        itCount = itCount+1;
    end
    if (Converged == 0)
        dt = dt/10;
        chops = chops + 1;
    end
end
%Stats
FIM.Iter(Ndt) = itCount-1;
FIM.Chops(Ndt) = chops;
Timers.Construct = TimerConstruct;
Timers.Solve = TimerSolve;
%Reshape S before quitting
S = reshape(s,Nx,Ny);

%Update production and injection data in wells
Inj.water(Ndt+1) = Inj.water(Ndt) + dt*Inj.PI*(K(1,Inj.x,Inj.y)...
    *K(2, Inj.x, Inj.y))^0.5*(Inj.p-Pw(Inj.x,Inj.y))*Inj.Mw;
Prod.water(Ndt+1) = Prod.water(Ndt) - dt*Prod.PI*(K(1,Prod.x,Prod.y)...
    *K(2, Prod.x, Prod.y))^0.5*(Prod.p-Pw(Prod.x,Prod.y))*Mw((Prod.y-1)*Nx+Prod.x);
Prod.oil(Ndt+1) =  Prod.oil(Ndt) - dt*Prod.PI*(K(1,Prod.x,Prod.y)...
    *K(2, Prod.x, Prod.y))^0.5*(Prod.p-P(Prod.x,Prod.y))*Mo((Prod.y-1)*Nx+Prod.x);
end