%Build FIM Residual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 4 April 2016
%Last modified: 22 April 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD Residual
% Builds Residual vector

%Input variables:
%   Status0: initial status
%   Status: new status
%   dt: timestep
%   Mnw: nonwetting phase saturation
%   Mw: wetting phase saturation
%   UpWindNw: nonwettin phase upwind operator
%   UpWindW: wetting phase upwind operator
%   Inj: injection wells
%   Prod: production wells
%   K: permeability


%Output variables:
%   Residual:  residual
%   Tnw: nonwetting phase transmissibility matrix
%   Tw: wetting phase transmissibility matrix

function [Residual, Tnw, Tw] = FIMResidual(Status0, Status, Grid, dt, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, K)

%Create local variables
N = Grid.N;
pv = Grid.por*Grid.Volume;
p_old = Status0.p;
s_old = Status0.s;
p = Status.p;
s = Status.s;
pc = Status.pc;

%Accumulation Term
Ap = speye(N)*0; %It's still incompressible
AS = speye(N)*pv/dt;

%Transmissibility matrix
Tnw = TransmissibilityMatrix (Grid, UpWindNw, Mnw);
Tw = TransmissibilityMatrix (Grid, UpWindW, Mw);

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
Residual = [Rw; Rnw];
end

%% Transmissibility
function T  = TransmissibilityMatrix(Grid, UpWind, M)
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Grid.N;

%%%Transmissibility matrix construction
Tx = zeros(Nx+1, Ny);
Ty = zeros(Nx, Ny+1);
%Apply upwind operator
Mupx = UpWind.x*M;
Mupy = UpWind.y*M;
Mupx = reshape(Mupx, Nx, Ny);
Mupy = reshape(Mupy, Nx, Ny);
Tx(2:Nx,:)= Grid.Tx(2:Nx,:).*Mupx(1:Nx-1,:);
Ty(:,2:Ny)= Grid.Ty(:,2:Ny).*Mupy(:,1:Ny-1);
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
for i=1:length(Inj)
    c = Inj(i).cells;
    switch (Inj(i).type)
        case('RateConstrained')
            qw(c) = Inj(i).q;
        case('PressureConstrained')
            qnw(c) = Inj(i).Mo * Inj(i).PI .* K(c).* (Inj(i).p - p(c));
            qw(c) = Inj(i).Mw * Inj(i).PI .* K(c) .* (Inj(i).p - p(c));
    end    
end
%Producers
for i=1:length(Prod)
    c = Prod(i).cells;
    switch (Prod(i).type)
        case('RateConstrained')
            qw(c) = Mw(c)./(Mw(c)+Mnw(c)).*Prod(i).q;
            qnw(c) = Prod(i).q - qw(c);
        case('PressureConstrained')
            qnw(c) =  Mnw(c).* Prod(i).PI .* K(c).* (Prod(i).p - p(c));
            qw(c) =   Mw(c).* Prod(i).PI .* K(c) .* (Prod(i).p - p(c));
    end
end
end