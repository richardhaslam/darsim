%Fractional Flow Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fw, dFw, DdFw] = ComputeFractionalFlow(s, Fluid)
%Compute Mobilities
[Mw, Mo] = Mobilities(s, Fluid);
Mt = Mw+Mo;   %total mobility
%Flux function and derivatives
fw = Mw./Mt;
df = Derivative(s,Fluid);
Ddf = Derivative_2nd(s, Fluid);

end