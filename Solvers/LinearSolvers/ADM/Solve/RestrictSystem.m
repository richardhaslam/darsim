%ADM - Restrict system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%Last modified: 2 May 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J_c] = RestrictSystem(J, R, Pp, Ps, N, maxLevel)
% Restrict system to DLGR Grid
Jop = J(1:N, 1:N);
Jos = J(1:N, N+1:2*N);
Jwp = J(N+1:2*N, 1:N);
Jws = J(N+1:2*N, N+1:2*N);
Jopc = R(1).matrix*Jop*Pp(1).matrix;
Josc = R(1).matrix*Jos*Ps(1).matrix;
Jwpc = R(1).matrix*Jwp*Pp(1).matrix;
Jwsc = R(1).matrix*Jws*Ps(1).matrix;
for i=2:maxLevel
    Jopc = R(i).matrix*Jopc*Pp(i).matrix;
    Josc = R(i).matrix*Josc*Ps(i).matrix;
    Jwpc = R(i).matrix*Jwpc*Pp(i).matrix;
    Jwsc = R(i).matrix*Jwsc*Ps(i).matrix;
end
J_c = [Jopc, Josc; Jwpc, Jwsc];
end