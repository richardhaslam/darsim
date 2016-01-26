function Fluxes = ComputeWellFluxes(Fluxes, Well, p, K)
for i=1:length(Well)
    a = Well(i).cells.cells;
    Fluxes(a) = Fluxes(a) + Well(i).PI.*K(a).*(Well(i).p-p(a));
end
end