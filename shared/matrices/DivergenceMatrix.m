%Divergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, q, U, Tx, Ty] = DivergenceMatrix(Grid, P, K, Tx, Ty, Inj, Prod)
%Builds Upwind Flux matrix
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Grid.N;

%Compute total fluxes ([m^3/s])
U.x = zeros(Nx+1,Ny,1);
U.y = zeros(Nx,Ny+1,1);
U.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)).*Tx(2:Nx,:);
U.y(:,2:Ny)  = (P(:,1:Ny-1)-P(:,2:Ny)).*Ty(:,2:Ny);

%Wells: fluxes [m^3/s]
Wells.Fluxes=zeros(Nx,Ny,1);
Wells.fw=zeros(Nx,Ny,1);
Wells.Fluxes(Inj.x,Inj.y)=Inj.PI*K(1,Inj.x,Inj.y)*(Inj.p-P(Inj.x,Inj.y));
Wells.Fluxes(Prod.x, Prod.y)=Prod.PI*K(1,Prod.x,Prod.y)*(Prod.p-P(Prod.x,Prod.y));

q=reshape(Wells.Fluxes, N,1);
A = SaturationMatrix(Grid, U, q);
end