% Compute phase saturations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 22 June 2016
%Last Modified: 22 June 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute phase saturations
% Updates the saturation of the 2 phases based on mass balance.
% 2 equations are used:
%   - z1*rho_t = x11*rho_1*S_1 + x12*rho_2*S_2
%   - rho_t = rho_1*S_1 + rho_2*S_2
% If the system is one-phase the saturation is set to 1 (vapor) or 0 (liquid)

function Status = ComputePhaseSaturation(Status, SinglePhase)
Status.s = Update_S(Status.x1, Status.z, Status.rho);
Status.s(SinglePhase.onlyvapor == 1) = 1;
Status.s(SinglePhase.onlyliquid == 1) = 0;
end

function [s] = Update_S(x,z,rho)
%Two phase, two component saturation updater    
    s = rho(:,2).*(x(:,2) - z(:,1))./(rho(:,1).*(z(:,1)...
        - x(:,1)) + rho(:,2).*(x(:,2) - z(:,1))); 
end
