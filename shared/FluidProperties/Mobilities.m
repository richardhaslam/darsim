%Phase mobilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mw, Mo, dMw, dMo]=Mobilities(s, Fluid)
%Phase Mobilities
S = (s-Fluid.swc)/(1-Fluid.swc-Fluid.sor); % Rescale saturations
if (strcmp(Fluid.RelPerm,'Quadratic')==1)
    Mw = S.^2 /Fluid.muw;        % Water mobility
    Mo = (1-S).^2/Fluid.muo;     % Oil mobility
    dMw = (1-Fluid.swc-Fluid.sor)^(-1)*2*S/Fluid.muw;
    dMo = -(1-Fluid.swc-Fluid.sor)^(-1)*2*(1-S)/Fluid.muo;
else
    Mw=S/Fluid.muw;              %Water mobility
    Mo=(1-S)/Fluid.muo;          %Oil mobility
    dMw = ones(length(s),1)*(1-Fluid.swc-Fluid.sor)^(-1)/Fluid.muw;
    dMo = -ones(length(s),1)*(1-Fluid.swc-Fluid.sor)^(-1)/Fluid.muo;
end
end