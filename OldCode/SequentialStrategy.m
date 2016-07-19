%SEQUENTIAL STRATEGY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%Last modified: 16 May 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Status, Inj, Prod, dT, Converged, Timers, ImplicitSolver] = SequentialStrategy(Status, K, Grid, Fluid, Inj, Prod, Sequential, Ndt, maxdT)
%SEQUENTIAL STRATEGY
Grid.CFL = Sequential.CFL;
MaxExtIter = Sequential.MaxExtIter;
Tol = Sequential.Tol;

%Initialise variables
s0 = Status.s; %previus timestep solution

%Timers inside timesteps
ptimer = zeros(MaxExtIter,1);
btimer = zeros(MaxExtIter,1);
stimer = zeros(MaxExtIter,1);

%External Newton loop
Iter = 1; %External iterations counter
Converged = 0;
s = s0;
while (Converged==0 && Iter <= MaxExtIter)
    disp(['Outer-loop iteration: ' num2str(Iter)]);
    tstart1 = tic;
    %1. Solve flow equation for pressure and compute fluxes
    disp('Pressure solver')
    [p, U, Pc, Wells, Inj, Prod] = PressureSolver(Grid, Inj, Prod, Fluid, S, K);
    ptimer(Iter) = toc(tstart1);
    
    %2. Check mass balance
    tstart2 = tic;
    [Balance, U] = check2D(U, Grid, Wells);
    btimer(Iter) = toc(tstart2);
    
    %3. Compute timestep-size based on CFL
    if (Iter==1)
        dT = timestepping(Fluid, Grid);
        dT = min(dT, maxdT);
    end
    
    %4. Solve transport equation given the total velocity field
    tstart4 = tic;
    sold = s; %Last converged solution
    if (Balance==1)
        q = reshape(Wells.Fluxes, Grid.N,1);
        if (Sequential.ImpSat==0)
            [s] = ExplicitTransport(Fluid, Grid, s0, U, q, dT);
            Converged = 1;
        else
            disp('Transport solver');
            Sequential.ImplicitSolver.timestep = [Sequential.ImplicitSolver.timestep, Ndt];
            [s, qnw, qw, Sequential.ImplicitSolver, dT, Tconverged] = ...
                ImplicitTransport(Fluid, Grid, s0, sold, U, q, Sequential.ImplicitSolver, dT, K);
            if Tconverged == 0
                disp('Transport solver did not converge')
                break
            end
            disp('----------------------------------')
        end
    else
        disp('Mass balance not respected!!');
        break
    end
    stimer(Iter) = toc(tstart4);
    
    %5. Compute DeltaS to check convergence
    Delta = (Fluid.s - sold);
    DeltaNorm = norm(Delta);
    if (DeltaNorm < Tol || MaxExtIter==1)
        Converged = 1;
        %%%%% Save converged solution %%%%
        Status.p = p;
        Status.pc = reshape(Pc, Grid.N, 1);
        Status.s = s;
    end
    Iter = Iter + 1;
end


%Compute Nwetting and wetting phase fluxes for production curves
for i=1:length(Prod)
    c = Prod(i).cells;
    Prod(i).qw = sum(qw(c));
    Prod(i).qnw = sum(qnw(c));
end

 if (Sequential.ImpSat==1)
    ImplicitSolver = Sequential.ImplicitSolver;
 else
     ImplicitSolver.timestep = 0;
     ImplicitSolver.Chops = 0;
     ImplicitSolver.Newtons = 0;
 end
Timers.ptimer = ptimer;
Timers.btimer = btimer;
Timers.stimer = stimer;
end