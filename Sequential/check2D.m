%Check mass balance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Balance, U]=check2D(U, Grid, Wells)
%Checks mass balance in all cells
Nx = Grid.Nx;
Ny = Grid.Ny;
Balance = 1;
maxUx = max(max(U.x));
maxUy = max(max(U.y));
maxU = max(maxUx, maxUy);
for i=1:Nx
    for j=1:Ny
        Accum = U.x(i,j) - U.x(i+1,j) + U.y(i,j) - U.y(i,j+1) + Wells.Fluxes(i,j);
        if (abs(Accum/maxU) > 10^(-5))
            Balance = 0;
        end
    end
end
 
end