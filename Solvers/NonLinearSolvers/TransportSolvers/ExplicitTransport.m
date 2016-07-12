%Explicit transport solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = ExplicitTransport(Fluid, Grid, S, U, q, dt)
%Explicit Transport Solver
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Nx*Ny;
por = Grid.por;
pv = por*Grid.Volume; 
s = reshape(S,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. Assemble system matrix
A = SaturationMatrix(Grid,U,q);      % system matrix
T = spdiags(dt/pv*ones(N,1),0,N,N);    % dt/pv * Cell Fluxes and producer
B = T*A;
injector = max(q,0).*dt/pv;                 % injection flux*dt/pv

% Update saturation
[Mw, Mo] = Mobilities(s, Fluid);
Mt = Mw+Mo; 
fw = Mw./Mt;
s = s +(B*fw + injector);                   
end