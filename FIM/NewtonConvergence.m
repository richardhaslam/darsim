%Check Convergence of Newton's loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last modified: 3 May 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Converged = NewtonConvergence(iter, Residual, Delta, Tol, N, ADM)
Converged = 0;
%Evaluate norms
if ADM.active == 1
   Residual_c = RestrictResidual(Residual, ADM.Rest, N, ADM.level);
   Norm1 =  norm(Residual_c, inf);
   Norm2 = norm(Delta, inf);
else
    Norm1 =  norm(Residual, inf);
    Norm2 = norm(Delta, inf);
end
disp(['Iter ' num2str(iter) '    ' num2str(Norm1), '    ', num2str(Norm2)]);
%Check convergence
if (Norm1 < Tol && Norm2 < Tol*1e3)
   Converged = 1;
end
end