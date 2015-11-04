%Assign neighbours
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Grid = AssignNeighbours(Grid)
Nx = Grid.Nx;
Ny = Grid.Ny;
i = 1;
for j=1:Ny
    g = i + (j-1)*Nx;
    if j~=1 && j~=Ny
    Grid.Neighbours(g).indexes = [g+1, g-Nx, g+Nx];
    elseif j == 1
        Grid.Neighbours(g).indexes = [g+1, g+Nx];
    elseif j == Ny
        Grid.Neighbours(g).indexes = [g+1, g-Nx];
    end
end
i = Nx;
for j=1:Ny
    g = i + (j-1)*Nx;
    if j~=1 && j~=Ny
    Grid.Neighbours(g).indexes = [g-1, g-Nx, g+Nx];
    elseif j == 1
        Grid.Neighbours(g).indexes = [g-1, g+Nx];
    elseif j == Ny
        Grid.Neighbours(g).indexes = [g-1, g-Nx];
    end
end
j=1;
for i=2:Nx-1
    g = i + (j-1)*Nx;
        Grid.Neighbours(g).indexes = [g-1, g+1, g+Nx];
end
j=Ny;
for i=2:Nx-1
    g = i + (j-1)*Nx;
        Grid.Neighbours(g).indexes = [g-1, g+1, g-Nx];
end

for i =2:Nx-1
    for j=2:Ny-1
        g = i + (j-1)*Nx;
        Grid.Neighbours(g).indexes = [g-1, g+1, g-Nx, g+Nx];
    end
end
end