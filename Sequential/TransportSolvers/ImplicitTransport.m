function [Snew, ImplicitSolver, dt]=ImplicitTransport(Fluid, Grid, S0, Sold, U, q, ImplicitSolver, dt)
%Implicit Transport Solver
Nx=Grid.Nx;
Ny=Grid.Ny;
N=Nx*Ny;
por=Grid.por;
pv=por*Grid.Volume;   %Void Volume in each cell
tol=ImplicitSolver.tol;
MaxIter=ImplicitSolver.maxNewton;
%%Create saturation vectors
sold=reshape(Sold, N, 1); %last converged saturation  
s0=reshape(S0, N, 1);    %saturation at previous timestep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
converged=0;
chops=0;
while (converged==0)   %If it does not converge the timestep is chopped
    %Initial guess for Newton loop
    snew=sold;
    %Residual at first iteration
    [Mw, Mo]=Mobilities(snew, Fluid);
    Mt=Mw+Mo;   %total mobility
    %Flux function and derivatives
    fw=Mw./Mt;
    df=Derivative(snew,Fluid);
    Ddf=Derivative_2nd(snew, Fluid);
    A=SaturationMatrix(Grid,U,q);      % system matrix
    Rafael = 0;
    Residual=max(q,0)+A*fw-pv/dt*(snew-s0);
    
    % Initialise objects
    Norm=1;
    dS=zeros(N,1);
    
    % Newton-Raphson loop
    Newton=1; %Newton's iterations counter
    dS_crit = zeros(N,1);
    while ((Norm>tol && Newton<=MaxIter) || (Newton==1))
        %Compute dS at nu+1
        D=spdiags(pv/dt*ones(N,1),0,N,N);
        B=D-A*spdiags(df,0,N,N);
        
        % FLUX CORRECTION - MATTEO
        if (ImplicitSolver.fluxfunction==1)
            dS_old = dS;
            dS=B\Residual;
            if (min(dS_crit .* dS) < 0)
                dS_crit = dS_old .* ((dS_crit .* dS) < 0)...
                    - dS_old .* ((dS_crit .* dS) > 0);
                C = B-A*spdiags(Ddf.*dS_crit,0,N,N);
                dS=C\Residual;
            end
        else
            dS=B\Residual;
        end
        
        %Update Saturation and remove unphysical values
        s_old=snew;
        snew = s_old+dS;
        snew = min(snew,1);
        snew = max(snew,0);
        
        % FLUX CORRECTION - PATRICK
        if (ImplicitSolver.fluxfunction==2)
            Ddf=Derivative_2nd(snew, Fluid);
            Ddf_old=Derivative_2nd(s_old, Fluid);
            snew=snew.*(Ddf.*Ddf_old>=0)+0.5*(snew+sold).*(Ddf.*Ddf_old<0);
            dS=snew-s_old;
        end
        
        %Compute norm of dS
        Norm=norm(dS, inf);
        
        %Compute Residual at nu
        [Mw, Mo]=Mobilities(snew, Fluid);   % mobilities at iteration nu
        Mt=Mw+Mo;   %total mobility
        fw=Mw./Mt;
        df=Derivative(snew,Fluid);
        Ddf=Derivative_2nd(snew, Fluid);
        Residual=max(q,0)+A*fw-pv/dt*(snew-s0);
        
        %Increase iteration counter
        Newton=Newton+1;
    end
    
    if (Norm<=tol) 
        converged=1;
    else
        chops=chops+1;
        dt=round(dt/2);
    end
end

%%%%Output objects
ImplicitSolver.Chops = [ImplicitSolver.Chops, chops];
ImplicitSolver.Newtons = [ImplicitSolver.Newtons, Newton];
%New saturation matrix
Snew=reshape(snew,Nx,Ny,1);
end