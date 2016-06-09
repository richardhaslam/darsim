%Time-step selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dt = timestepping(Fluid, Grid, U)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns the timestep size based on the velocity field, the Grid and the
%fluid present. CFL condition is used. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
por = Grid.por;
pv = por*Grid.Volume;   %Void Volume in each cell
CFL = Grid.CFL;

%I take the worst possible scenario
s = Fluid.sr(2):0.01:1-Fluid.sr(1);
df =  ComputeFractionalFlow(s, Fluid);
dfmax = max(df);
Uxmax = max(max(abs(U.x)));
Uymax = max(max(abs(U.y)));
Lambdax = dfmax * Uxmax;
Lambday = dfmax * Uymax;

%Compute timestep size
dtx = CFL*pv/Lambdax;
dty = CFL*pv/Lambday;
dt = min(dtx,dty);


end