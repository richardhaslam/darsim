%Fluid Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fluid = FluidProperties(density, viscosity, relperm, compressibility, capillarity, foam, Comp_Type, Comp_Prop, inputMatrix)
%Fluid Properties
% Composition properties
Fluid.Type=char(inputMatrix(Comp_Type + 1));
Fluid.NumPhase=str2double(inputMatrix(Comp_Type + 3));
Fluid.Comp.NumComp=str2double(inputMatrix(Comp_Type + 5));
% RelPerm function: Quadratic or Linear
Fluid.RelPerm = char(inputMatrix(relperm + 1));
% Capillary pressure function: polynomial/exponential/table
Fluid.Pc = char(inputMatrix(capillarity + 1));  
fplot = 0; %hard-coded for now

    for getDataPhase = 1:Fluid.NumPhase
       %Gets all densities [kg/m^3]
       Fluid.rho(getDataPhase)=str2double(inputMatrix(density + 2*getDataPhase));
       %Gets all viscosities [Pa sec]
       Fluid.mu(getDataPhase)=str2double(inputMatrix(viscosity + 2*getDataPhase));
       %Gets all Corey exponents [-]
       Fluid.n(getDataPhase)=str2double(inputMatrix(relperm + 9 + 2*getDataPhase));
       %Gets all compressibilities [1/Pa]
       Fluid.c(getDataPhase)=str2double(inputMatrix(compressibility + 2*getDataPhase));
       %Gets all residual saturations [-]
       Fluid.sr(getDataPhase)=str2double(inputMatrix(relperm + 1 + 2*getDataPhase));
       %Gets all end-point rel perms [-]
       Fluid.kre(getDataPhase)=str2double(inputMatrix(relperm + 5 + 2*getDataPhase));
    end
if isnan(Fluid.Comp.NumComp)
else
    for getDataComp=1:Fluid.Comp.NumComp
       %Gets all atmospheric bubble points [K]
       Comp_Prop + 3 + (getDataComp-1)*5
       Fluid.Comp.Tb(getDataComp)=str2double(inputMatrix(Comp_Prop + 3 + (getDataComp-1)*5));
       %Gets all slopes connecting bubble point and critical point on 1/T plot
       %[K]
       Fluid.Comp.b(getDataComp)=str2double(inputMatrix(Comp_Prop + 5*getDataComp));
    end
end

% Foam properties
Fluid.epdry=str2double(inputMatrix(foam + 2));
Fluid.fmdry=str2double(inputMatrix(foam + 4));
Fluid.fmmob=str2double(inputMatrix(foam + 6));
    
if (fplot==1)
    FractionalFlowPlot(Fluid);
end
end