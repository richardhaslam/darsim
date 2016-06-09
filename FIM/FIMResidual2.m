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
%   p_old: previous timestep pressure solution
%   s_old: previous timestep saturation solution
%   p: new pressure
%   s: new saturaton
%   pc: capillary pressure
%   pv: pore volume
%   dt: timestep
%   Trx: rock transmissibility in x direction
%   Try: rock transmissibility in y direction
%   Mnw: nonwetting phase saturation
%   Mw: wetting phase saturation
%   UpWindNw: nonwettin phase upwind operator
%   UpWindW: wetting phase upwind operator
%   Inj: injection wells
%   Prod: production wells
%   K: permeability
%   N: number of grid cells
%   Nx: number of grid cells in x direction
%   Ny: number of grid cells in y direction

%Output variables:
%   Residual:  residual
%   Tnw: nonwetting phase transmissibility matrix
%   Tw: wetting phase transmissibility matrix

function [Residual, qtot, Tnw, Tw] = FIMResidual2(Status0, Status, pc, pv, dt, Trx, Try, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod, K, Grid, Fluid)
p_old = Status0.p;
s_old = Status0.s;
x1_old = Status0.x1;
p = Status.p;
s = Status.s;
x1 = Status.x1;
N = Grid.N;
Nx = Grid.Nx;
Ny = Grid.Ny;

%Density
[Rho,~] = LinearDensity(p,Fluid.c,Fluid.rho);
[Rho_old,~] = LinearDensity(p_old,Fluid.c,Fluid.rho);

%Accumulation Term
Ap = speye(N)*0; %It's still incompressible
AS = speye(N)*pv/dt;

%Transmissibility matrix
Tnw = TransmissibilityMatrix (Trx, Try, N, Nx, Ny, UpWindNw, Mnw);
Tw = TransmissibilityMatrix (Trx, Try, N, Nx, Ny, UpWindW, Mw);

%Gravity
G = ComputeGravityTerm(N);

%Source terms
[qnw, qw] = ComputeWells(N, Inj, Prod, K, Status, Mnw, Mw);
qtot = qnw + qw;

%% RESIDUAL
%Non-wetting phase
Rnw = Ap*p - Ap*p_old...                                                %Compressibility term
    - AS*s.*x1(:,1).*Rho(:,1) + AS*s_old.*x1_old(:,1).*Rho_old(:,1)...                      %Accumulation due to phase 1
    + AS*s.*x1(:,2).*Rho(:,2) - AS*s_old.*x1_old(:,2).*Rho_old(:,2)...                      %Accumulation due to phase 2
    + Tw*(p-pc).*x1(:,2).*Rho(:,2) + Tnw*p.*x1(:,1).*Rho(:,1)...                            %Convection
    + G*s - qnw.*x1(:,1).*Rho(:,1) - qw.*x1(:,2).*Rho(:,2);                                 %Gravity and source terms
%Wetting phase (Has capillarity in immiscible model)
Rw = Ap*p - Ap*p_old...                                                 %Compressibility term
    - AS*s.*(1-x1(:,1)).*Rho(:,1) + AS*s_old.*(1-x1_old(:,1)).*Rho_old(:,1)...              %Accumulation due to phase 1
    + AS*s.*(1-x1(:,2)).*Rho(:,2) - AS*s_old.*(1-x1_old(:,2)).*Rho_old(:,2)...              %Accumulation due to phase 2
    + Tw*(p-pc).*(1-x1(:,2)).*Rho(:,2) + Tnw*p.*(1-x1(:,1)).*Rho(:,1)...                    %Convection
    + G*s - qnw.*(1-x1(:,1)).*Rho(:,1) - qw.*(1-x1(:,2)).*Rho(:,2);                         %Gravity and source terms
%Stick them together
Residual = [Rnw; Rw];
Residual = Residual/Rho(1,1);
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
function [qnw, qw] = ComputeWells(N, Inj, Prod, K, Status, Mnw, Mw)
p = Status.p;
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