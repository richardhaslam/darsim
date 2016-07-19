% Compute z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 22 June 2016
%Last Modified: 22 June 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute z

function [z] = Compute_Z(s, x, rho)
%Two phase, two component total mole fraction updater
z(:,1) = (x(:,1).*rho(:,1).*s + x(:,2).*...
    rho(:,2).*(1-s))./(rho(:,1).*s + rho(:,2).*(1-s));  %Based on mass balance equation z_1 * rho_t = x11*rho1*s1 + x12*rho2*s2
z(:,2) = 1 - z(:,1);

end

