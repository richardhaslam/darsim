%Construct upwind operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last modified: 5 April 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upwind operator
% constructs upwind operator for a given phase pressure

%Input variables:
%   Grid: grid information
%   P: phase pressure 

%Output variables:
%   A: upwind operator
%   U: velocity

function [A, U] = UpwindOperator (Grid, P)
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Grid.N;

%Compute 'rock' fluxes ([m^3/s])
U.x = zeros(Nx+1,Ny,1);
U.y = zeros(Nx,Ny+1,1);
U.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)).*Grid.Tx(2:Nx,:);
U.y(:,2:Ny)  = (P(:,1:Ny-1)-P(:,2:Ny)).*Grid.Ty(:,2:Ny);

%Use velocity to build upwind operator
R = reshape((U.x(2:Nx+1,:) >= 0), N, 1);
L = reshape((U.x(1:Nx,:) < 0), N, 1);
T = reshape((U.y(:,2:Ny+1) >= 0), N, 1);
B = reshape((U.y(:,1:Ny) < 0), N, 1);
DiagVecs = [R, L];
DiagIndx = [0,1];
A.x = spdiags(DiagVecs,DiagIndx,N,N);
DiagVecs = [T, B];
DiagIndx = [0,Nx];
A.y = spdiags(DiagVecs,DiagIndx,N,N);
end