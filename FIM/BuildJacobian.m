%Build FIM Jacobian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last Modified: 5 April 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD Jacobian
% Builds Jacobian matrix

%Input variables:
%   Grid: grid
%   K: permeability 
%   TMatrixNw: nonwetting phase transmissibility matrix
%   TMatrixW: wetting phase transmissibility matrix
%   p: pressure soltuion at previous iteration
%   Mw: wetting phase mobility 
%   Mnw: nonwetting phase mobility 
%   dMw: derivative of wetting phase mobility 
%   dMnw: derivative of nonwetting phase mobility  
%   Unw: non-wetting phase rock flux 
%   Uw: wetting phase rock flux
%   dPc: derivative of capillary pressure 
%   dt: timestep size
%   Inj: injection wells 
%   Prod: production wells
%   UpWindNw: upwind operator of nonwetting phase
%   UpWindW: upwind operator of wetting phase

%Output variables:
%   J:  FIM Jacobian

function J = BuildJacobian(Grid, K, TMatrixNw, TMatrixW, p, Mw, Mnw, dMw, dMnw, Unw, Uw, dPc, dt, Inj, Prod, UpWindNw, UpWindW)
%Build FIM Jacobian
Nx = Grid.Nx; 
Ny = Grid.Ny; 
N = Grid.N;
pv = Grid.Volume*Grid.por;

% BUILD FIM JACOBIAN BLOCK BY BLOCK

%1. Rnw Pressure Block
Jnwp = TMatrixNw;

%2. Rw Pressure Block
Jwp = TMatrixW;

%3. Rnw Saturation Block
dMupxNw = UpWindNw.x*dMnw;
dMupyNw = UpWindNw.y*dMnw;
%Construct Jnws block
x1 = min(reshape(Unw.x(1:Nx,:),N,1),0).*dMupxNw;
x2 = max(reshape(Unw.x(2:Nx+1,:),N,1),0).*dMupxNw;
y1 = min(reshape(Unw.y(:,1:Ny),N,1),0).*dMupyNw;
y2 = max(reshape(Unw.y(:,2:Ny+1),N,1),0).*dMupyNw;
v = ones(N,1)*pv/dt;
DiagVecs = [-y2, -x2, y2+x2-y1-x1-v, x1, y1];
DiagIndx = [-Nx, -1, 0, 1, Nx];
JnwS = spdiags(DiagVecs,DiagIndx,N,N);

%4. Rw Saturation Block
dMupxw = UpWindW.x*dMw;
dMupyw = UpWindW.y*dMw;
%Construct JwS block
x1 = min(reshape(Uw.x(1:Nx,:),N,1),0).*dMupxw;
x2 = max(reshape(Uw.x(2:Nx+1,:),N,1),0).*dMupxw;
y1 = min(reshape(Uw.y(:,1:Ny),N,1),0).*dMupyw;
y2 = max(reshape(Uw.y(:,2:Ny+1),N,1),0).*dMupyw;
v = ones(N,1)*pv/dt;
DiagVecs = [-y2, -x2, y2+x2-y1-x1+v, x1, y1];
DiagIndx = [-Nx, -1, 0, 1, Nx];
JwS = spdiags(DiagVecs,DiagIndx,N,N);
CapJwS = full(Jwp);

%Add capillarity
for i = 1:N
     CapJwS(i,:) = CapJwS(i,:).* dPc';
end
JwS = JwS - sparse(CapJwS);

% Add wells
[Jnwp, Jwp, JnwS, JwS] = AddWellsToJacobian(Jnwp, Jwp, JnwS, JwS, Inj, Prod, K, p, Mw, Mnw, dMw, dMnw);

% Full Jacobian: put the 4 blocks together
J = [Jnwp, JnwS; Jwp, JwS];
end