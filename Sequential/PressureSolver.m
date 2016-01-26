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
Kvector = reshape(K(1,:,:), N, 1);

%Transmissibility
[Tx, Ty] = ComputeTransmissibility(Grid, K);

%Construct pressure matrix
A = AssemblePressureMatrix(Tx, Ty, Nx, Ny);
Ab = A;
q=zeros(N,1);

%Add Wells
for i=1:length(Inj)
    [A, q] = AddWell(A, q, Inj(i), Kvector);
end
for i=1:length(Prod)
    [A, q] = AddWell(A, q, Prod(i), Kvector);
end
%Solve for pressure
p = A\q;

%Compute total fluxes ([m^3/s])
P = reshape(p, Nx, Ny,1);
U.x = zeros(Nx+1,Ny,1);
U.y = zeros(Nx,Ny+1,1);
U.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)).*Tx(2:Nx,:);
U.y(:,2:Ny)  = (P(:,1:Ny-1)-P(:,2:Ny)).*Ty(:,2:Ny);

%Wells: fluxes [m^3/s]
Fluxes = zeros(N,1);
Fluxes = ComputeWellFluxes(Fluxes, Inj, p, Kvector); 
Fluxes = ComputeWellFluxes(Fluxes, Prod, p, Kvector);
Wells.Fluxes = reshape(Fluxes, Nx, Ny);
end