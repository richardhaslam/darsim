function [Rho, dRhodP] = LinearDensity( P, cf, Rho0)
%Calculates the density of a fluid at a new pressure when given the
%original density at a reference pressure. Does the same for porosity when
%required.

Prho = 1e5;
Rho = zeros(length(P),length(cf));
dRhodP = zeros(length(P),length(cf));
%needs to be vectorized
%Rho = Rho0.*exp(cf.*(P - Prho));
%dRhodP = cf.*Rho0.*exp(cf.*(P(i)-Prho));
for i = 1:length(P)  
    Rho(i,:) = Rho0.*exp(cf.*(P(i)-Prho));
    dRhodP(i,:) = cf.*Rho0.*exp(cf.*(P(i)-Prho));
end

end

