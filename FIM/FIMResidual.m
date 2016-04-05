%Build FIM Residual for one phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 4 April 2016
%Last modified: 5 April 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Residual, Tnw, Tw] = FIMResidual(p_old, s_old, p, s, pc, pv, dt, Trx, Try, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, K, N, Nx, Ny)
%Accumulation Term
Ap = speye(N)*0; %It's still incompressible
AS = speye(N)*pv/dt;

%Transmissibility matrix
Tnw = TransmissibilityMatrix (Trx, Try, N, Nx, Ny, UpWindNw, Mnw);
Tw = TransmissibilityMatrix (Trx, Try, N, Nx, Ny, UpWindW, Mw);

%Gravity

%Source terms
q = zeros(N,1);
q = ComputeWellFluxes(q, Inj, p, K, p*0, K); 
q = ComputeWellFluxes(q, Prod, p, K, p*0, K);
qnw = Inj(1).Mo * max(q,0) + Mnw.*min(q, 0);
qw = Inj(1).Mw * max(q,0) + Mw.*min(q, 0);

%% RESIDUAL
%Non-wetting phase
Rnw = Ap*(p-p_old) - AS*(s-s_old) + Tnw*p - qnw;
%Wetting phase (This one has capillarity)
Rw = Ap*(p-p_old) + AS*(s-s_old) + Tw*p - Tw*pc - qw;
%Stick them together
Residual = [Rnw; Rw];

end


function T  = TransmissibilityMatrix(Trx, Try, N, Nx, Ny, UpWind, M)
%%%Transmissibility matrix construction
Tx = zeros(Nx+1, Ny);
Ty = zeros(Nx, Ny+1);
%Apply upwind operator
Mupx = UpWind.x*M;
Mupy = UpWind.y*M;
Mupx = reshape(Mupx, Nx, Ny);
Mupy = reshape(Mupy, Nx, Ny);
Tx(2:Nx,:)= Trx(2:Nx,:).*Mupx(1:Nx-1,:);
Ty(:,2:Ny)= Try(:,2:Ny).*Mupy(:,1:Ny-1);
%Construct matrix 
x1 = reshape(Tx(1:Nx,:),N,1);
x2 = reshape(Tx(2:Nx+1,:),N,1);
y1 = reshape(Ty(:,1:Ny),N,1);
y2 = reshape(Ty(:,2:Ny+1),N,1);
DiagVecs = [-y2,-x2,y2+x2+y1+x1,-x1,-y1];
DiagIndx = [-Nx,-1,0,1,Nx];
T = spdiags(DiagVecs,DiagIndx,N,N);
end