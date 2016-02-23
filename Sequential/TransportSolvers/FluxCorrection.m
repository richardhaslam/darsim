%Flux Correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [snew, dS] = FluxCorrection(snew, sold, Fluid)
% FLUX CORRECTION - PATRICK
if (ImplicitSolver.fluxfunction==2)
    Ddf = Derivative_2nd(snew, Fluid);
    Ddf_old = Derivative_2nd(s_old, Fluid);
    snew = snew.*(Ddf.*Ddf_old>=0)+0.5*(snew+sold).*(Ddf.*Ddf_old<0);
    dS = snew-s_old;
else 
    dS = snew - sold;
end
end