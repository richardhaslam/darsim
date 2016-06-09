%Check Input File Errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Grid.Tres < 273.15
    display('Are you sure that reservoir temperature is in Kelvin?')
    Errors = 1;
end
if Grid.por < 0 || Grid.por > 1
    display('Invalid porosity input')
    Errors = 1;
end
if max(Fluid.rho) < 100
    display('Most likely wrong units for fluid density (must be kg/m^3)')
    Errors = 1;
end
if min(Fluid.c) < 0
    display('Invalid fluid compressibility')
    Errors = 1;
end
if max(Fluid.kre) > 1
    display('Invalid end-point relative permeability')
    Errors = 1;
end
if isempty(Fluid.Type) == 1
    display('Please enter compositional model type (Immiscible,BlackOil,Compositional)')
    Errors = 1;
end
if (strcmp(Fluid.Type,'Compositional')==1)
    if isempty(Fluid.Comp.Tb) == 1 || max(isnan(Fluid.Comp.Tb)) == 1
        display('Please enter the atmospheric boiling point of the components [K]')
        Errors = 1;
    end
    if isempty(Fluid.Comp.b) == 1 || max(isnan(Fluid.Comp.b)) == 1
        display('Please enter the b parameter of the components [K]')
        Errors = 1;
    end
end
if (strcmp(Fluid.Type,'Immiscible')==1) || (strcmp(Fluid.Type,'BlackOil')==1) || (strcmp(Fluid.Type,'Compositional')==1)
else
    display('Please enter a valid compositional model (Immiscible,BlackOil,Compositional)')
    Errors = 1;
end
if (strcmp(Fluid.Type,'BlackOil')==1) || (strcmp(Fluid.Type,'Compositional')==1)
    if isempty(FlashSettings.TolInner) == 1 || max(isnan(FlashSettings.TolInner)) == 1
        display('Please enter a convergence tolerance for the inner loop (flash) update')
        Errors = 1;
    end
    if isempty(FlashSettings.MaxIt) == 1 || max(isnan(FlashSettings.MaxIt)) == 1
        display('Please enter a maximum number of iterations for the inner loop (flash) update')
        Errors = 1;
    end
    if(strcmp(Fluid.Type,'Compositional')==1)
        if isempty(FlashSettings.TolFlash) == 1 || max(isnan(FlashSettings.TolFlash)) == 1
            display('Please enter a tolerance for the compositional flash update')
            Errors = 1;
        end
    end
end
if isempty(Fluid.NumPhase) == 1
    display('Please enter the number of phases')
    Errors = 1;
end
if(strcmp(Fluid.Type,'Compositional')==1)
    if isempty(Fluid.Comp) == 1
        display('Please enter the number of components')
        Errors = 1;
    end
    if isempty(Grid.Tres) == 1
        display('Please enter a reservoir temperature')
        Errors = 1;
    end
end
if isempty(Grid.por) == 1
    display('Please enter a porosity')
    Errors = 1;
end
if isempty(Fluid.RelPerm) == 1
    display('Please enter the type of rel perm curve to be used (Quadratic,Foam,Other)')
    Errors = 1;
end
if (strcmp(Fluid.RelPerm,'Quadratic')==1) || (strcmp(Fluid.RelPerm,'Foam')==1) || (strcmp(Fluid.RelPerm,'Linear')==1) || (strcmp(Fluid.RelPerm,'Corey')==1)
else
    display('Please enter a valid rel perm model (Quadratic,Foam,Linear,Corey)')
    Errors = 1;
end
if (strcmp(Fluid.Pc,'Linear')==1) || (strcmp(Fluid.Pc,'JLeverett')==1) || (strcmp(Fluid.Pc,'BrooksCorey')==1) || (strcmp(Fluid.Pc,'Table')==1)
elseif isempty(Fluid.Pc) == 1
else
    display('Please enter capillarity model type (Linear,JLeverett,BrooksCorey,Table) or delete any reference to a capillary model')
    Errors = 1;
end

if isempty(Fluid.rho) == 1 || max(isnan(Fluid.rho)) == 1
    display('Please enter fluid densities [kg/m^3]')
    Errors = 1;
end
if isempty(Fluid.mu) == 1 || max(isnan(Fluid.mu)) == 1
    display('Please enter viscosities [Pa sec]')
    Errors = 1;
end
if isempty(Fluid.c) == 1 || max(isnan(Fluid.c)) == 1
    display('Please enter fluid compressibilities [1/Pa]')
    Errors = 1;
end
if isempty(Fluid.sr) == 1 || max(isnan(Fluid.sr)) == 1
    display('Please enter residual saturations [-]')
    Errors = 1;
end
if (strcmp(Fluid.RelPerm,'Foam')==1) || (strcmp(Fluid.RelPerm,'Corey')==1)
    if isempty(Fluid.kre) == 1 || max(isnan(Fluid.kre)) == 1
        display('Please enter end-point relative permeabilities [-]')
        Errors = 1;
    end
    if isempty(Fluid.n) == 1 || max(isnan(Fluid.n)) == 1
        display('Please enter Corey function exponents [-]')
        Errors = 1;
    end
end
if (strcmp(Fluid.RelPerm,'Foam')==1)
    if isempty(Fluid.epdry) == 1 || isempty(Fluid.fmdry) == 1 || isempty(Fluid.fmmob) == 1
        display('Please enter foam properties (epdry,fmdry,fmmob)')
        Errors = 1;
    end
end