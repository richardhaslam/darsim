%ADM - Restrict residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last modified: 21 March 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Residual_c = RestrictResidual(Residual, R, N, maxLevel)
Residual_cp = R(1).matrix*Residual(1:N);
Residual_cs = R(1).matrix*Residual(N+1:2*N);
for i=2:maxLevel
    Residual_cp = R(i).matrix*Residual_cp;
    Residual_cs = R(i).matrix*Residual_cs;
end
Residual_c = [Residual_cp; Residual_cs];
end