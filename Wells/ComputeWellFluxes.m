%Well Fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fluxes = ComputeWellFluxes(Fluxes, Well, p, K, pc, Kw)
for i=1:length(Well)
    a = Well(i).cells.cells;
    Fluxes(a) = Fluxes(a) + Well(i).PI.* (K(a).*(Well(i).p-p(a)) + Kw(a).*pc(a));
end
end