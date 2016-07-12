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

function Converged = NewtonConvergence(iter, Residual, Delta, p, Tol, N, ADM, Delta_c)
Converged = 0;
%Evaluate norms
if ADM.active == 1
    Residual_c = RestrictResidual(Residual, ADM.Rest, N, ADM.level);
    Nc = length(Delta_c)/2;
    Norm1 =  norm(Residual_c, inf);
    Norm2 = norm(Delta_c(1:Nc), inf)/max(p);
    Norm3 = norm(Delta_c(Nc+1:end), inf);
else
    Norm1 =  norm(Residual, inf);
    Norm2 = norm(Delta(1:N), inf)/max(p);
    Norm3 = norm(Delta(N+1:end), inf);
end

disp(['Iter ' num2str(iter) '    ' num2str(Norm1, '%5.5e'), '    ', num2str(Norm2,'%5.5e'), '    ', num2str(Norm3, '%5.5e')]);

%Check convergence
if (Norm1 < Tol && Norm2<Tol && Norm3 < Tol)
    Converged = 1;
end
end