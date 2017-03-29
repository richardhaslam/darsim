%FIM Linear Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last modified: 21 March 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Linear Solver
%Given a linear system A x = b it solves it either with ADM or without

%Input variables:
%   A: matrix
%   RHS: right-hand side 
%   N: number of fine-scale gridblocks
%   ADM: contains ADM objects

%Output variables:
%   x: solution

function [x, x_c] = LinearSolver (A, RHS, N, ADM)
if ADM.active == 1 && ADM.level~=0
    RHS_c = RestrictResidual(RHS, ADM.Rest, N, ADM.level);
    [A_c] = RestrictSystem(A, ADM.Rest, ADM.Prolp, ADM.Prols, N, ADM.level);
    x_c = -A_c\RHS_c;
    x = Prolong(x_c, ADM.Prolp, ADM.Prols, ADM.level);
else
    x = -A\RHS;
    x_c = x;
end
end