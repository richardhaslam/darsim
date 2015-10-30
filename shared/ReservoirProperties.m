function [Grid, Inj, Prod, Lx, Ly, K] = ReservoirProperties(Problem)
%Define Reservoir Properties, Grid and wells
%Dimensions
Lx = 540;                              %Dimension in x−direction [m] 
Ly = 540;                              %Dimension in y−direction [m]
h = 1;                                  %Reservoir thickness [m]

%Gridding
Grid.Nx = 54; 
Grid.dx = Lx/Grid.Nx; 
Grid.Ny = 54; 
Grid.dy = Ly/Grid.Ny; 
Grid.Ax = Grid.dy*h;                    %Cross section in x direction
Grid.Ay = Grid.dx*h;                    %Cross section in y direction
Grid.Volume = Grid.dx.*Grid.dy*h;       %Cell volume [m^3]
Grid.por = 0.2;                         %Porosity

%Rock permeability in [m^2].
K = ReadPerm(Grid, Problem);

%%%%Wells%%%%
%Injection wells
Inj.r = 0.15; %Well radius in m
Inj.p = 10^5; %[Pa]
Inj.x = 1;
Inj.y = Grid.Ny;
Inj.PI = ComputeProductivityIndex(Inj.r, K(1, Inj.x, Inj.y), K(2, Inj.x, Inj.y), Grid.dx, Grid.dy, 1);
%Production Wells
Prod.r = 0.15; %Well radius in m
Prod.p = 0; %[Pa]
Prod.x = Grid.Nx;
Prod.y = 1;
Prod.PI = ComputeProductivityIndex(Prod.r, K(1, Prod.x, Prod.y), K(2, Prod.x, Prod.y), Grid.dx, Grid.dy, 1);
%%Plot Permeability Field
PlotPermeability(K, Grid, Lx, Ly);
end