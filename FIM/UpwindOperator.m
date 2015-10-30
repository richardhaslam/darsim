%Construct upwind operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A] = UpwindOperator (Grid, U)
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Nx*Ny;
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