%Compute Capillary Pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pc, dPc] = ComputePc(S, Fluid)
%Compute capillary pressure
%   Given the wetting phase saturation it returns Pc and dPc
switch(Fluid.Pc)
    case('polynomial')
        %Define parameters
        a=300;
        b=10;
        c=350;
        d=40;
        e=100;
        %Compute Pc and dPc analytically
        Pc = a + b*(1-S) + c*(1-S).^2 + d*(1-S).^3 +e*(1-S).^4;
        dPc = - b - 2*c*(1-S) - d*3*(1-S).^2 - e*4*(1-S).^3;
    case('BrooksCorey')
        %Define parameters
        Pct = 700;
        lambda = 2.5;
        %Compute Pc and dPc analytically
        Pc = Pct.*((1-Fluid.swc)./(S-Fluid.swc)).^(1/lambda);
        dPc = Pct*((1 - Fluid.swc)./(S - Fluid.swc)).^(1/lambda)./(lambda*(Fluid.swc - S));
    case('BentsenAnli')
        %Define parameters
        Pct = 700;
        Pcs = 350;
        %Compute Pc and dPc analytically
        Pc = Pct - Pcs.*log((S-Fluid.swc)/(1-Fluid.swc));
        dPc = Pcs./(Fluid.swc - S);
    case('Table')
        %here you read the table
end
end

