function Fluxes = ComputeWellFluxes(Fluxes, Well, P, K)
n = length(Well.x);
for i=1:n
    Fluxes(Well.x(i),Well.y(i)) = Well.PI*K(1,Well.x(i),Well.y(i))*(Well.p-P(Well.x(i),Well.y(i)));
end
end