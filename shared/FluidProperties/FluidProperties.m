%Fluid Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fluid = FluidProperties(viscosity, relperm, capillarity, foam, corey, inputMatrix)
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
% Endpoint relperms
Fluid.krwe=str2double(inputMatrix(corey + 2));
Fluid.kroe=str2double(inputMatrix(corey + 4));
% Corey exponents
Fluid.nw=str2double(inputMatrix(corey + 6));
Fluid.no=str2double(inputMatrix(corey + 8));
% Foam properties
Fluid.epdry=str2double(inputMatrix(foam + 2));
Fluid.fmdry=str2double(inputMatrix(foam + 4));
Fluid.fmmob=str2double(inputMatrix(foam + 6));
fplot = 0; %hard-coded for now
if (fplot==1)
    FractionalFlowPlot(Fluid);
end
if ~isempty(Fluid.Pc)
    PcPlot(Fluid);
end
end