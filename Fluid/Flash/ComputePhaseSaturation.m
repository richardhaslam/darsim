% Compute phase saturations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 22 June 2016
%Last Modified: 22 June 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute pahse saturations
% Updates the saturation of the 2 phases based on mass balance.
% 2 options:
%   a. black oil: 
%                  
%   b. Compositional: 


function Status = ComputePhaseSaturation(Status, SinglePhase, rho)
Status.s = Update_S(Status.x1, Status.z, rho);
Status.s(SinglePhase.onlyvapor == 1) = 1;
Status.s(SinglePhase.onlyliquid == 1) = 0;
end

function [s] = Update_S(x,z,rho)
%Two phase, two component saturation updater    
    s = rho(:,2).*(x(:,2) - z(:,1))./(rho(:,1).*(z(:,1)...
        - x(:,1)) + rho(:,2).*(x(:,2) - z(:,1)));  %Based on mass balance equation found in Update_z.m  
end
