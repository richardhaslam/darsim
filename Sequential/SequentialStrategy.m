%SEQUENTIAL STRATEGY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, S, Pc, dT, Converged, Timers, ImplicitSolver] = SequentialStrategy(S0, K, Grid, Fluid, Inj, Prod, Sequential, Ndt, maxdT)
%SEQUENTIAL STRATEGY
Grid.CFL = Sequential.CFL;
MaxExtIter = Sequential.MaxExtIter;
N = Grid.Nx*Grid.Ny;
Tol = Sequential.Tol;

%Timers inside timesteps
ptimer = zeros(MaxExtIter,1);
btimer = zeros(MaxExtIter,1);
stimer = zeros(MaxExtIter,1);
%External Newton loop
Iter = 1; %External iterations counter
Converged = 0;
S = S0;
while (Converged==0 && Iter <= MaxExtIter)
    disp(['Outer-loop iteration: ' num2str(Iter)]);
    tstart1 = tic;
    %1. Solve flow equation for pressure and compute fluxes
    disp('Pressure solver')
    [P, U, Pc, Wells] = PressureSolver(Grid, Inj, Prod, Fluid, S, K);
    %% 
    ptimer(Iter) = toc(tstart1);
    
    %2. Check mass balance
    tstart2 = tic;
    [Balance, U] = check2D(U, Grid, Wells);
    btimer(Iter) = toc(tstart2);
    
    %3. Compute timestep-size based on CFL
    if (Iter==1)
        dT = timestepping(Fluid, S, Grid, U, Wells);
        dT = min(dT, maxdT);
        dT = 1500/30*24*3600;
    end
    
    %4. Solve transport equation given the total velocity field
    tstart4 = tic;
    Sold=S; %Last converged solution
    if (Balance==1)
        q=reshape(Wells.Fluxes, N,1);
        if (Sequential.ImpSat==0)
            [S]=ExplicitTransport(Fluid, Grid, S0, U, q, dT);
            Converged = 1;
        else
            disp('Transport solver');
            Sequential.ImplicitSolver.timestep = [Sequential.ImplicitSolver.timestep, Ndt];
            [S, Sequential.ImplicitSolver, dT, Tconverged] = ImplicitTransport(Fluid, Grid, S0, Sold, U, q, Sequential.ImplicitSolver, dT, K);
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
    Delta=(S-Sold);
    DeltaNorm = norm(reshape(Delta, N,1));
    if (DeltaNorm < Tol || MaxExtIter==1)
        Converged = 1;
    end
    Iter = Iter+1;
end
 if (Sequential.ImpSat==1)
    ImplicitSolver=Sequential.ImplicitSolver;
 else
     ImplicitSolver.timestep = 0;
     ImplicitSolver.Chops = 0;
     ImplicitSolver.Newtons = 0;
 end
Timers.ptimer = ptimer;
Timers.btimer = btimer;
Timers.stimer = stimer;
end
