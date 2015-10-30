%Check mass balance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Balance=check2D(U, Grid, Wells)
%Checks mass balance in all cells
Nx=Grid.Nx;
Ny=Grid.Ny;
Balance=1;
for i=1:Nx
    for j=1:Ny
        Accum=U.x(i,j)-U.x(i+1,j)+U.y(i,j)-U.y(i,j+1)+Wells.Fluxes(i,j);
        if (Accum>10^(-10))
            Balance=0;
        end
    end
end
 
end