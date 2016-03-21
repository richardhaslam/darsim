%FIM Linear Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last modified: 21 March 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Delta] = LinearSolver (J, Residual, N, ADM)
if ADM.active == 1
    Residual_c = RestrictResidual(Residual, ADM.Rest, N, ADM.level);
    [J_c] = RestrictSystem(J, ADM.Rest, ADM.Prolp, ADM.Prols, N, ADM.level);
    Delta_c = -J_c\Residual_c;
    Delta = Prolong(Delta_c, ADM.Prolp, ADM.Prols, ADM.level);
else
    Delta = -J\Residual;
end
end