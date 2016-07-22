%Flux Correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [snew, dS] = FluxCorrection(snew, sold, Fluid, fluxcorrection)
% FLUX CORRECTION - PATRICK
if (fluxcorrection == 2)
    Ddf = Derivative_2nd(snew, Fluid);
    Ddf_old = Derivative_2nd(sold, Fluid);
    snew = snew.*(Ddf.*Ddf_old>=0)+0.5*(snew+sold).*(Ddf.*Ddf_old<0);
    dS = snew-sold;
else 
    dS = snew - sold;
end
end