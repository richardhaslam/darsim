%Fractional Flow Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fw, dfw] = ComputeFractionalFlow(s, Fluid)
%Compute Mobilities
[Mw, Mnw, dMw, dMnw] = Mobilities(s, Fluid);
Mt = Mw + Mnw;
%Flux function and derivatives
fw = Mw./Mt;
dfw = (dMw .* (Mw+Mnw) - (dMw+dMnw).* Mw)./(Mw+Mnw).^2;
end