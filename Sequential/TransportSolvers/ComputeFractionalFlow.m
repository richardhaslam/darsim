%Fractional Flow Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fw, dfw] = ComputeFractionalFlow(s, Fluid)
%Compute Mobilities
[Mw, Mo] = Mobilities(s, Fluid);
Mt = Mw + Mo;
%Flux function and derivatives
fw = Mw./Mt;
dfw = Derivative(s,Fluid);
end