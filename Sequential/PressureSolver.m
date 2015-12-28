%Pressure Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, U, Wells, A, Ab, q]=PressureSolver(Grid, K, Inj, Prod, Mt)
%PRESSURE Solver

%1.Compute transmissibilities using harmonic average.
Nx = Grid.Nx; Ny = Grid.Ny; 
N = Grid.N;

%Transmissibility
[Tx, Ty] = ComputeTransmissibility(Grid, K);
% Tx (2:Nx+1,:) = Tx(2:Nx+1,:).*Mt;
% Ty (:,2:Ny+1) = Ty(:,2:Ny+1).*Mt;

%Construct pressure matrix
A = AssemblePressureMatrix(Tx, Ty, Nx, Ny);
Ab = A;
q=zeros(N,1);

%Add Wells
[A, q] = AddWell(A, q, Inj, K, Nx);
[A, q] = AddWell(A, q, Prod, K, Nx);

%Solve for pressure
p = A\q;

%Compute total fluxes ([m^3/s])
P = reshape(p, Nx, Ny,1);
U.x = zeros(Nx+1,Ny,1);
U.y = zeros(Nx,Ny+1,1);
U.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)).*Tx(2:Nx,:);
U.y(:,2:Ny)  = (P(:,1:Ny-1)-P(:,2:Ny)).*Ty(:,2:Ny);

%Wells: fluxes [m^3/s]
Wells.Fluxes=zeros(Nx,Ny);
Wells.Fluxes = ComputeWellFluxes(Wells.Fluxes, Inj, P, K); 
Wells.Fluxes = ComputeWellFluxes(Wells.Fluxes, Prod, P, K);
end