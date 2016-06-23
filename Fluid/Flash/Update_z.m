
function [z] = Update_z(s, x, rho)
%Two phase, two component total mole fraction updater
z(:,1) = (x(:,1).*rho(:,1).*s(:,1) + x(:,2).*...
    rho(:,2).*s(:,2))./(rho(:,1).*s(:,1) + rho(:,2).*s(:,2));  %Based on mass balance equation z_1 * rho_t = x11*rho1*s1 + x12*rho2*s2
z(:,2) = 1 - z(:,1);

end

