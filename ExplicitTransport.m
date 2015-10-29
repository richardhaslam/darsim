function [Snew]=ExplicitTransport(Fluid, Grid, S, U, q, dt)
%Explicit Transport Solver
Nx=Grid.Nx;
Ny=Grid.Ny;
N=Nx*Ny;
por=Grid.por;
pv=por*Grid.Volume;   %Void Volume in each cell
s=reshape(S,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. Assemble system matrix
A=SaturationMatrix(Grid,U,q);      % system matrix
T=spdiags(dt/pv*ones(N,1),0,N,N);    % dt/pv * Cell Fluxes and producer
B=T*A;
fi=max(q,0).*dt/pv;                 % injection flux*dt/pv

% Update saturation
[Mw, Mo]=Mobilities(s, Fluid);
Mt=Mw+Mo;   %total mobility
fw=Mw./Mt;
s = s +(B*fw+fi);                   
Snew=reshape(s,Nx,Ny,1);
end
