%Add Pc to right-hand side of the pressure equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q, Pc, U] = AddPcToPressureSystem(q, S, Fluid, Kw, K, Grid)
N = Grid.Nx*Grid.Ny;
Nx = Grid.Nx;
Ny = Grid.Ny;
Q = reshape(q, Nx, Ny);
[Tx, Ty] = ComputeTransmissibility(Grid, Kw);

%Compute Pc for the Saturation distribution given
[Pc, ~] = ComputePc(S, Fluid, K, Grid.por);

%Add Pc to the right-hand side
%"Capillary fluxes"
U.x = zeros(Nx+1,Ny,1);
U.y = zeros(Nx,Ny+1,1);
U.x(2:Nx,:) = (Pc(1:Nx-1,:)-Pc(2:Nx,:)).*Tx(2:Nx,:);
U.y(:,2:Ny)  = (Pc(:,1:Ny-1)-Pc(:,2:Ny)).*Ty(:,2:Ny);
%Right-hand side
Q (:,:) = Q(:,:) - U.x(1:Nx, :) + U.x(2:Nx+1,:) - U.y(:,1:Ny) + U.y(:,2:Ny+1);
%reshape
q = reshape(Q, N, 1);
end