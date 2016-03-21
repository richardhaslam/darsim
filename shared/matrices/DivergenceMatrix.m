%Divergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, Fluxes, U] = DivergenceMatrix(Grid, P, K, Tx, Ty, Inj, Prod)
%Builds Upwind Flux matrix
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Grid.N;
Kvector = reshape(K(1,:,:), N, 1);

%Compute total fluxes ([m^3/s])
U.x = zeros(Nx+1,Ny,1);
U.y = zeros(Nx,Ny+1,1);
U.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)).*Tx(2:Nx,:);
U.y(:,2:Ny)  = (P(:,1:Ny-1)-P(:,2:Ny)).*Ty(:,2:Ny);

%Wells: fluxes [m^3/s]
Fluxes = zeros(N,1);
p = reshape(P, N, 1);
Fluxes = ComputeWellFluxes(Fluxes, Inj, p, Kvector, p*0, Kvector); 
Fluxes = ComputeWellFluxes(Fluxes, Prod, p, Kvector, p*0, Kvector);
Wells.Fluxes = reshape(Fluxes, Nx, Ny);

A = SaturationMatrix(Grid, U, Fluxes);
end