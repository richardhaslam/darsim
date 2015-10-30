%MsFV solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p_f = MsFVSolve(A, q, G, MsR, MsP)
% MS pressure solution
% Permute
q = G*q;
A = G*A*G';
% Restrict
A_c = MsR*A*MsP;
q_c = MsR*(q - A*C*q);

% Solve
p_c = A_c\q_c;

% Prolong
p_f = G'*(MsP*p_c) + G'*C*G*q;
end