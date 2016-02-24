%Implicit transport solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Snew, ImplicitSolver, dt, converged]=ImplicitTransport(Fluid, Grid, S0, Sold, U, q, ImplicitSolver, dt, K)
%Implicit Transport Solver
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Nx*Ny;
por = Grid.por;
pv = por*Grid.Volume;   %Void Volume in each cell
tol = ImplicitSolver.tol;
MaxIter = ImplicitSolver.maxNewton;
%%Create saturation vectors
sold = reshape(Sold, N, 1); %last converged saturation  
s0 = reshape(S0, N, 1);    %saturation at previous timestep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
converged = 0;
chops = 0;
while (converged==0 && chops<=10)   %If it does not converge the timestep is chopped
    %Initial guess for Newton loop
    snew = sold;
    
    %Residual at first iteration
    help =  K(1, :, :);
    K = zeros(Nx, Ny);
    K(:,:) = help(1, :,:);
    Residual = TransportResidual(snew, s0, q, pv, U, dt, Fluid, Grid, K);
    
    % Initialise objects
    Norm = 1;
    Newton = 1; %Newton's iterations counter
    
    % Newton-Raphson loop
    while ((Norm > tol && Newton <= MaxIter) || (Newton==1))
        %Compute dS at nu+1
        D = spdiags(pv/dt*ones(N,1),0,N,N);
        B = D - A * spdiags(df,0,N,N);
        dS = B \ Residual;
        
        %Update Saturation and remove unphysical values
        sold = snew;
        snew = sold + dS;
        snew = min(snew,1);
        snew = max(snew,Fluid.swc);
        
        [snew, dS] = FluxCorrection(snew, sold, Fluid, ImplicitSolver.fluxfunction);
        
        %Compute norm of dS
        Norm = norm(dS, inf);
        
        %Compute Residual at nu
        TransportResidual(snew, s0, q, pv, U, dt, Fluid, Grid)
        
        %Increase iteration counter
        Newton=Newton+1;
    end
    if (Norm <= tol) 
        converged = 1;
    else
        chops = chops + 1;
        dt = dt/2;
    end
end

%%%%Output objects
ImplicitSolver.Chops = [ImplicitSolver.Chops, chops];
ImplicitSolver.Newtons = [ImplicitSolver.Newtons, Newton];
%New saturation matrix
Snew=reshape(snew,Nx,Ny,1);

end