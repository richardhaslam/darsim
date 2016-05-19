%Check Convergence of Newton's loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last modified: 3 May 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Newton convergence check
% Checks if the non-linear convergence has been obtained

%Input variables:
%   iter: iteration number
%   Residual:   residual 
%   Delta: correction to the solution
%   Tol: tolerance
%   N: number of fine-scale gridblocks
%   ADM: contains ADM objects

%Output variables:
%   Converged:  1 if it's a converged solution, 0 otherwise

function Converged = NewtonConvergence(iter, Residual, qtot, Delta, Tol, N, ADM)
Converged = 0;
%Evaluate norms
if ADM.active == 1
    Residual_c = RestrictResidual(Residual, ADM.Rest, N, ADM.level);
    Norm1 =  norm(Residual_c, inf) / max(abs(qtot));
    Norm2 = norm(Delta, inf);
else
    Norm1 =  norm(Residual, inf) / max(abs(qtot));
    Norm2 = norm(Delta, inf);
end
disp(['Iter ' num2str(iter) '    ' num2str(Norm1), '    ', num2str(Norm2)]);
%Check convergence
if (Norm1 < Tol && Norm2 < Tol*1e3)
   Converged = 1;
end
end