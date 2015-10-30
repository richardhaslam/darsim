%Fluid Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fluid = FluidProperties(fplot)
%Fluid Properties
Fluid.muw=1*10^(-5);  
Fluid.muo=1.5*10^(-5); % Viscosities [Pa*s]
Fluid.swc=0; Fluid.sor=0;           % Irreducible saturations
Fluid.RelPerm = 'Quadratic';            % RelPerm function: Quadratic or Linear
if (fplot==1)
    FractionalFlowPlot(Fluid);
end
end