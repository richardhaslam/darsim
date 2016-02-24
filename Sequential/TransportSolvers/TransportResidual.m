%Compute Transport Residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Residual, A, df] = TransportResidual(snew, s0, q, pv, U, dt, Fluid, Grid, K)
%Compute residual
[fw, df] = ComputeFractionalFlow(snew, Fluid);
A = SaturationMatrix(Grid,U,q);      % Total fluxes matrix
B = CapFluxMatrix(Grid, reshape(snew, Grid.Nx, Grid.Ny), Fluid, K);
Residual = pv/dt*(snew-s0)  - max(q,0) - A*fw - B*ones(Grid.N,1);
end

function B = CapFluxMatrix(Grid, S, Fluid, K)
%Builds a matrix that sums capillary fluxes
CapFlux.x = zeros(Grid.Nx+1, Grid.Ny);
CapFlux.y = zeros(Grid.Ny, Grid.Ny +1);
C = ComputeCancesFunction(Fluid,S, K, Grid);
CapFlux.x(2:Nx, :) = C(1:Nx-1, :) - C(2:Nx, :);
CapFlux.y(:,2:Ny) = C(:,1:Ny-1) - C(:, 2:Ny);
%
x1 = reshape(CapFlux.x(1:Nx,:),N,1);  
y1 = reshape(CapFlux.y(:,1:Ny),N,1);
% − flow in positive coordinate (XP, YP, ZP)
x2 = reshape(CapFlux.x(2:Nx+1,:),N,1);
y2 = reshape(CapFlux.y(:,2:Ny+1),N,1);

%assemble matrix
DiagVecs = -x1 - y1 + x2 + y2; % diagonal vectors
DiagIndx = 0; % diagonal index
B = spdiags(DiagVecs,DiagIndx,N,N);
end