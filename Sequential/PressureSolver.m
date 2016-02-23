%Pressure Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, U, Wells, A, Ab, q]=PressureSolver(Grid, Inj, Prod, Fluid, S, K)
%PRESSURE Solver

%1.Compute transmissibilities using harmonic average.
Nx = Grid.Nx; Ny = Grid.Ny; 
N = Grid.N;

%Transmissibility
% Effective permeability
[Mw, Mo]=Mobilities(S, Fluid);
Mt=Mw+Mo;   %total mobility
Kt=zeros(2, Grid.Nx, Grid.Ny);
Kt(1,:,:)=reshape(Mt, 1, Grid.Nx, Grid.Ny).*K(1,:,:);		% x-direction
Kt(2,:,:)=reshape(Mt, 1, Grid.Nx, Grid.Ny).*K(2,:,:);		% y-direction
Kvector = reshape(Kt(1,:,:), N, 1);
[Tx, Ty] = ComputeTransmissibility(Grid, Kt);

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

%Add capillary term to the right-hand side
if ~isempty(Fluid.Pc)
    Kw(1,:,:)=reshape(Mw, 1, Grid.Nx, Grid.Ny).*K(1,:,:);		% x-direction
    Kw(2,:,:)=reshape(Mw, 1, Grid.Nx, Grid.Ny).*K(2,:,:);		% y-direction
    [q, Pc] = AddPcToPressureSystem(q, S, Fluid, Kw, Grid);
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