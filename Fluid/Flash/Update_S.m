function [s] = Update_S(x,z,rho)
%Two phase, two component saturation updater    
    s(:,1) = rho(:,2).*(x(:,2) - z(:,1))./(rho(:,1).*(z(:,1)...
        - x(:,1)) + rho(:,2).*(x(:,2) - z(:,1)));  %Based on mass balance equation found in Update_z.m
    s(:,2) = 1 - s(:,1);    
end

