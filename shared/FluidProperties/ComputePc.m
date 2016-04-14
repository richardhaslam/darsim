%Compute Capillary Pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pc, dPc, dJ] = ComputePc(S, Fluid, K, por)
%Compute capillary pressure
%   Given the wetting phase saturation it returns Pc and dPc
switch(Fluid.Pc)
    case ('Linear')
        S = (S - Fluid.swc)./(1 - Fluid.swc);
        PcScale = 1e5; %in Pa
        Pc = PcScale.*(1-S);
        dPc = -PcScale;
        dJ = dPc;
    case('JLeverett')
        %J-leverett curve
        S = (S - Fluid.swc)./(1 - Fluid.swc);
        J = 0.1*(S).^(-0.5);
        dJ = - 0.1*0.5*(S).^(-1.5);
        %Define parameters
        alpha = 4.361e-5;
        %Compute Pc and dPc analytically
        Pc = alpha * (por./K).^(0.5).*J;
        dPc = alpha * (por./K).^(0.5).* dJ;
    case('BrooksCorey')
        %Define parameters
        Pct = 700;
        lambda = 2.5;
        %Compute Pc and dPc analytically
        Pc = Pct.*((1-Fluid.swc)./(S-Fluid.swc)).^(1/lambda);
        dPc = Pct*((1 - Fluid.swc)./(S - Fluid.swc)).^(1/lambda)./(lambda*(Fluid.swc - S));
        dJ = dPc;
    case('Table')
        %here you read the table
        disp('ERROR: The option <Table> for capillary pressure has not been implemented yet')
        return
    otherwise
        Pc = 0.*S;
        dPc = Pc;
end
end

