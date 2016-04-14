%Build FIM Residual 
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
G = ComputeGravityTerm(N);

%Source terms
[qnw, qw] = ComputeWells(N, Inj, Prod, K, p, Mnw, Mw);

%% RESIDUAL
%Non-wetting phase
Rnw = Ap*(p-p_old) - AS*(s-s_old) + Tnw*p + G*s - qnw;
%Wetting phase (This one has capillarity)
Rw = Ap*(p-p_old) + AS*(s-s_old) + Tw*p - Tw*pc + G*s - qw;
%Stick them together
Residual = [Rnw; Rw];
end

%% Transmissibility
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

%% Gravity
function G = ComputeGravityTerm(N)
G = speye(N)*0;
end
%% Wells
function [qnw, qw] = ComputeWells(N, Inj, Prod, K, p, Mnw, Mw)
qnw = zeros(N,1);
qw = zeros(N,1);
%Injectors
for i=1:size(Inj)
    c = Inj.cells;
    qnw(c) = Inj(i).Mo * Inj(i).PI .* K(c).* (Inj(i).p - p(c));
    qw(c) = Inj(i).Mw * Inj(i).PI .* K(c) .* (Inj(i).p - p(c));
end
%Producers
for i=1:size(Prod)
    c = Prod(i).cells;
    qnw(c) =  Mnw(c).* Prod(i).PI .* K(c).* (Prod(i).p - p(c));
    qw(c) =   Mw(c).* Prod(i).PI .* K(c) .* (Prod(i).p - p(c));
end
qt = qnw + qw;
qnw(1);
qnw(end);
qw(1);
qw(end);
qt(1);
qt(end);
end