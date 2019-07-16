%*********************************************************************
%*                                                                   *
%*                    Micro Boundary constraints                     *
%*                                                                   *
%*********************************************************************
%
% This function is usefull to impose the boundary conditions on the micro
% scale problem.
% Remark: The boundary condition will be Periodic or dirichlet, the user
% can handle it.
% Case Dirichlet: We impose weakle the Dirichlet boundary condition like an
% new data of the liner system.
% Periodic: We impose the peridicity of the soluton at the boundary using
% the selected pairs of edges (Pre-process) and we impose the condition
% like a constraints of the problem, for this reason we use 'lsqlin' to
% solve linear systems with constraints.

function [Sol1,Sol2] = Micro_Boundary(M,b1,b2,Micro_geo)

%% Auxiliar matrix for the periodicity

neq = Micro_geo.noedges+Micro_geo.nElement;

Indic_periodic = sparse([1:neq Micro_geo.edgePar(:,2)'],...
    [1:neq Micro_geo.edgePar(:,1)'],[ones(neq,1);-ones(size(Micro_geo.edgePar,1),1)]);
Indic_periodic(:,Micro_geo.edgePar(:,2))=[];

%% ZERO AVERAGE
last = [zeros(Micro_geo.noedges,1); Micro_geo.area'];
M  = [M last; last' 0];
b1 = [b1; 0]; b2 = [b2; 0];

Indic_periodic = [Indic_periodic; [zeros(1,size(Indic_periodic,2)-1),1]];
Indic_periodic = [Indic_periodic [zeros(size(Indic_periodic,1)-1,1);1]];

MM = (Indic_periodic'*M*Indic_periodic);
b1 = (Indic_periodic'*b1);
b2 = (Indic_periodic'*b2);

% -----------------------------------------------------------
%            SOLUTION OF DE EDGE PROBLEM
% -----------------------------------------------------------

x1 = MM\b1;
x2 = MM\b2;

% Copy the solution
indSolut = setdiff(1:neq+1,Micro_geo.edgePar(:,2));
Sol1 = zeros(neq+1,1);
Sol1(indSolut,1) = x1;
Sol1(Micro_geo.edgePar(:,2),1) =  -Sol1(Micro_geo.edgePar(:,1),1);

Sol2 = zeros(neq+1,1);
Sol2(indSolut,1) = x2;
Sol2(Micro_geo.edgePar(:,2),1) =  -Sol2(Micro_geo.edgePar(:,1),1);

end




