%ADM - Prolong solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Delta = Prolong(Delta_c, Pp, Ps, maxLevel)
%Prolongs pressure Solution from DLGR grid to fine Grid
Nc = length(Delta_c)/2;
Deltap = Pp(maxLevel).matrix*Delta_c(1:Nc);
Deltas = Ps(maxLevel).matrix*Delta_c(Nc+1:2*Nc);
for i = maxLevel-1:-1:1
    Deltap = Pp(i).matrix*Deltap;
    Deltas = Ps(i).matrix*Deltas;
end
Delta = [Deltap; Deltas];
end