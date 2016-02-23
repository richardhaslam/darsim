%Compute Transport Residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Residual = TransportResidual(snew, s0, q, Grid, U)
A = SaturationMatrix(Grid,U,q);      % system matrix
Residual = max(q,0) + A*fw - pv/dt*(snew-s0);     
end