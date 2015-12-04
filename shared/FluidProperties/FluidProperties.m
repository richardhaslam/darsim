%Fluid Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fluid = FluidProperties(viscosity, relperm, capillarity, inputMatrix)
%Fluid Properties
% Viscosities [Pa*s]
Fluid.muw=str2double(inputMatrix(viscosity + 2));  
Fluid.muo=str2double(inputMatrix(viscosity + 4)); 
% RelPerm function: Quadratic or Linear
Fluid.RelPerm = char(inputMatrix(relperm + 1));
% Capillary pressure function: polynomial/exponential/table
Fluid.Pc = char(inputMatrix(capillarity + 1));  
% Irreducible saturations
Fluid.swc=str2double(inputMatrix(relperm + 3));  
Fluid.sor=str2double(inputMatrix(relperm + 5));
fplot = 0; %hard-coded for now
if (fplot==1)
    FractionalFlowPlot(Fluid);
end
end