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

function J = BuildJacobian2(Grid, K, TMatrixNw, TMatrixW, Status, Mw, Mnw, dMw, dMnw, Unw, Uw, dPc, dt, Inj, Prod, UpWindNw, UpWindW, Fluid, pv)
p = Status.p;
s = Status.s;
x1 = Status.x1;
[Rho,dRhodP] = LinearDensity(p,Fluid.c,Fluid.rho);

%Build FIM Jacobian
Nx = Grid.Nx; 
Ny = Grid.Ny; 
N = Grid.N;
pv = Grid.Volume*Grid.por;

%% BUILD FIM JACOBIAN BLOCK BY BLOCK

%% 1. Saturation derivative in convection term
%For component 1
dMupxNw = UpWindNw.x*(dMnw.*Rho(:,1).*x1(:,1) + dMw.*Rho(:,2).*x1(:,2));
dMupyNw = UpWindNw.y*(dMnw.*Rho(:,1).*x1(:,1) + dMw.*Rho(:,2).*x1(:,2));
%Construct J1S block
left = min(reshape(Unw.x(1:Nx,:),N,1),0).*dMupxNw;
right = max(reshape(Unw.x(2:Nx+1,:),N,1),0).*dMupxNw;
bottom = min(reshape(Unw.y(:,1:Ny),N,1),0).*dMupyNw;
top = max(reshape(Unw.y(:,2:Ny+1),N,1),0).*dMupyNw;
DiagVecs = [-top, -right, top+right-bottom-left, left, bottom];
DiagIndx = [-Nx, -1, 0, 1, Nx];
J1S = spdiags(DiagVecs,DiagIndx,N,N);

%For component 2
dMupxW = UpWindW.x*(dMnw.*Rho(:,1).*(1-x1(:,1)) + dMw.*Rho(:,2).*(1-x1(:,2)));
dMupyW = UpWindW.y*(dMnw.*Rho(:,1).*(1-x1(:,1)) + dMw.*Rho(:,2).*(1-x1(:,2)));
%Construct J2S block
left = min(reshape(Uw.x(1:Nx,:),N,1),0).*dMupxW;
right = max(reshape(Uw.x(2:Nx+1,:),N,1),0).*dMupxW;
bottom = min(reshape(Uw.y(:,1:Ny),N,1),0).*dMupyW;
top = max(reshape(Uw.y(:,2:Ny+1),N,1),0).*dMupyW;
DiagVecs = [-top, -right, top+right-bottom-left, left, bottom];
DiagIndx = [-Nx, -1, 0, 1, Nx];
J2S = spdiags(DiagVecs,DiagIndx,N,N);


%% 2. Pressure derivative in convection term due to compressibility
%For component 1
dRhopxNw = UpWindNw.x*(Mnw.*dRhodP(:,1).*x1(:,1)) + UpWindW.x*(Mw.*dRhodP(:,2).*x1(:,2));
dRhopyNw = UpWindNw.y*(Mnw.*dRhodP(:,1).*x1(:,1)) + UpWindW.y*(Mw.*dRhodP(:,2).*x1(:,2));
%Construct J1p block
left = min(reshape(Unw.x(1:Nx,:),N,1),0).*dRhopxNw;
right = max(reshape(Unw.x(2:Nx+1,:),N,1),0).*dRhopxNw;
bottom = min(reshape(Unw.y(:,1:Ny),N,1),0).*dRhopyNw;
top = max(reshape(Unw.y(:,2:Ny+1),N,1),0).*dRhopyNw;
DiagVecs = [-top, -right, top+right-bottom-left, left, bottom];
DiagIndx = [-Nx, -1, 0, 1, Nx];
J1P = spdiags(DiagVecs,DiagIndx,N,N);

%For component 2
dRhopxW = UpWindNw.x*(Mnw.*dRhodP(:,1).*(1-x1(:,1))) + UpWindW.x*(Mw.*dRhodP(:,2).*(1-x1(:,2)));
dRhopyW = UpWindNw.y*(Mnw.*dRhodP(:,1).*(1-x1(:,1))) + UpWindW.y*(Mw.*dRhodP(:,2).*(1-x1(:,2)));
%Construct J2p block
left = min(reshape(Uw.x(1:Nx,:),N,1),0).*dRhopxW;
right = max(reshape(Uw.x(2:Nx+1,:),N,1),0).*dRhopxW;
bottom = min(reshape(Uw.y(:,1:Ny),N,1),0).*dRhopyW;
top = max(reshape(Uw.y(:,2:Ny+1),N,1),0).*dRhopyW;
DiagVecs = [-top, -right, top+right+bottom+left, -left, -bottom];
DiagIndx = [-Nx, -1, 0, 1, Nx];
J2P = spdiags(DiagVecs,DiagIndx,N,N);


%% 3. Pressure derivative in convection term due to flow
%For component 1 due to nw phase
dPpxNw = UpWindNw.x*(Rho(:,1).*x1(:,1));
dPpyNw = UpWindNw.y*(Rho(:,1).*x1(:,1));
%Construct J1P block dependence on nw phase
left(reshape(Unw.x(1:Nx,:),N,1)>=0,1) = dPpxNw(reshape(Unw.x(1:Nx,:),N,1)>=0,1);      %Flow from left block influence
left(reshape(Unw.x(1:Nx,:),N,1)<0,1) = 0;                                             %Flow flowing into left block so it has no influence
right(reshape(Unw.x(2:Nx+1,:),N,1)>0,1) = 0;                                           %Flow flowing into right block so it has no influence
right(reshape(Unw.x(2:Nx+1,:),N,1)<=0,1) = dPpxNw(reshape(Unw.x(2:Nx+1,:),N,1)<=0,1);  %Flow from right block influence
bottom(reshape(Unw.y(:,1:Ny),N,1)>=0,1) = dPpyNw(reshape(Unw.y(:,1:Ny),N,1)>=0,1);      %Flow from bottom block influence
bottom(reshape(Unw.y(:,1:Ny),N,1)<0,1) = 0;                                             %Flow flowing into bottom block so it has no influence
top(reshape(Unw.y(:,2:Ny+1),N,1)>0,1) = 0;                                           %Flow flowing into top block so it has no influence
top(reshape(Unw.y(:,2:Ny+1),N,1)<=0,1) = dPpyNw(reshape(Unw.y(:,2:Ny+1),N,1)<=0,1);  %Flow from top block influence
DiagVecs = [-top, -left, zeros(N,1), -right, -bottom];
DiagIndx = [-Nx, -1, 0, 1, Nx];
J1P = J1P - spdiags(DiagVecs,DiagIndx,N,N).*TMatrixNw;

%For component 1 due to w phase
dPpxW = UpWindW.x*(Rho(:,2).*x1(:,2));
dPpyW = UpWindW.y*(Rho(:,2).*x1(:,2));
%Construct J1P block dependence on w phase
left(reshape(Uw.x(1:Nx,:),N,1)>=0,1) = dPpxW(reshape(Uw.x(1:Nx,:),N,1)>=0,1);      %Flow from left block influence
left(reshape(Uw.x(1:Nx,:),N,1)<0,1) = 0;                                             %Flow flowing into left block so it has no influence
right(reshape(Uw.x(2:Nx+1,:),N,1)>0,1) = 0;                                           %Flow flowing into right block so it has no influence
right(reshape(Uw.x(2:Nx+1,:),N,1)<=0,1) = dPpxW(reshape(Uw.x(2:Nx+1,:),N,1)<=0,1);  %Flow from right block influence
bottom(reshape(Uw.y(:,1:Ny),N,1)>=0,1) = dPpyW(reshape(Uw.y(:,1:Ny),N,1)>=0,1);      %Flow from bottom block influence
bottom(reshape(Uw.y(:,1:Ny),N,1)<0,1) = 0;                                             %Flow flowing into bottom block so it has no influence
top(reshape(Uw.y(:,2:Ny+1),N,1)>0,1) = 0;                                           %Flow flowing into top block so it has no influence
top(reshape(Uw.y(:,2:Ny+1),N,1)<=0,1) = dPpyW(reshape(Uw.y(:,2:Ny+1),N,1)<=0,1);  %Flow from top block influence
DiagVecs = [-top, -left, zeros(N,1), -right, -bottom];
DiagIndx = [-Nx, -1, 0, 1, Nx];
J1P = J1P - spdiags(DiagVecs,DiagIndx,N,N).*TMatrixW;
dummy = sum(J1P);
J1P = J1P - diag(dummy);

%For component 2 due to nw phase
dPpxNw = UpWindNw.x*(Rho(:,1).*(1-x1(:,1)));
dPpyNw = UpWindNw.y*(Rho(:,1).*(1-x1(:,1)));
%Construct J2P block dependence on nw phase
left(reshape(Unw.x(1:Nx,:),N,1)>=0,1) = dPpxNw(reshape(Unw.x(1:Nx,:),N,1)>=0,1);      %Flow from left block influence
left(reshape(Unw.x(1:Nx,:),N,1)<0,1) = 0;                                             %Flow flowing into left block so it has no influence
right(reshape(Unw.x(2:Nx+1,:),N,1)>0,1) = 0;                                           %Flow flowing into right block so it has no influence
right(reshape(Unw.x(2:Nx+1,:),N,1)<=0,1) = dPpxNw(reshape(Unw.x(2:Nx+1,:),N,1)<=0,1);  %Flow from right block influence
bottom(reshape(Unw.y(:,1:Ny),N,1)>=0,1) = dPpyNw(reshape(Unw.y(:,1:Ny),N,1)>=0,1);      %Flow from bottom block influence
bottom(reshape(Unw.y(:,1:Ny),N,1)<0,1) = 0;                                             %Flow flowing into bottom block so it has no influence
top(reshape(Unw.y(:,2:Ny+1),N,1)>0,1) = 0;                                           %Flow flowing into top block so it has no influence
top(reshape(Unw.y(:,2:Ny+1),N,1)<=0,1) = dPpyNw(reshape(Unw.y(:,2:Ny+1),N,1)<=0,1);  %Flow from top block influence
DiagVecs = [-top, -left, zeros(N,1), -right, -bottom];
DiagIndx = [-Nx, -1, 0, 1, Nx];
J2P = J2P - spdiags(DiagVecs,DiagIndx,N,N).*TMatrixNw;

%For component 2 due to w phase
dPpxW = UpWindW.x*(Rho(:,2).*(1-x1(:,2)));
dPpyW = UpWindW.y*(Rho(:,2).*(1-x1(:,2)));
%Construct J2P block dependence on w phase
left(reshape(Uw.x(1:Nx,:),N,1)>=0,1) = dPpxW(reshape(Uw.x(1:Nx,:),N,1)>=0,1);      %Flow from left block influence
left(reshape(Uw.x(1:Nx,:),N,1)<0,1) = 0;                                             %Flow flowing into left block so it has no influence
right(reshape(Uw.x(2:Nx+1,:),N,1)>0,1) = 0;                                           %Flow flowing into right block so it has no influence
right(reshape(Uw.x(2:Nx+1,:),N,1)<=0,1) = dPpxW(reshape(Uw.x(2:Nx+1,:),N,1)<=0,1);  %Flow from right block influence
bottom(reshape(Uw.y(:,1:Ny),N,1)>=0,1) = dPpyW(reshape(Uw.y(:,1:Ny),N,1)>=0,1);      %Flow from bottom block influence
bottom(reshape(Uw.y(:,1:Ny),N,1)<0,1) = 0;                                             %Flow flowing into bottom block so it has no influence
top(reshape(Uw.y(:,2:Ny+1),N,1)>0,1) = 0;                                           %Flow flowing into top block so it has no influence
top(reshape(Uw.y(:,2:Ny+1),N,1)<=0,1) = dPpyW(reshape(Uw.y(:,2:Ny+1),N,1)<=0,1);  %Flow from top block influence
DiagVecs = [-top, -left, zeros(N,1), -right, -bottom];
DiagIndx = [-Nx, -1, 0, 1, Nx];
J2P = J2P - spdiags(DiagVecs,DiagIndx,N,N).*TMatrixW;
dummy = sum(J2P);
J2P = J2P - diag(dummy);

%% 4. Saturation derivative in accumulation term
%For component 1
v = ones(N,1)*pv/dt.*(x1(:,1).*Rho(:,1) - x1(:,2).*Rho(:,2));
J1S = J1S + spdiags(-v,0,N,N);
%For component 2
v = ones(N,1)*pv/dt.*((1-x1(:,1)).*Rho(:,1) - (1-x1(:,2)).*Rho(:,2));
J2S = J2S + spdiags(-v,0,N,N);


%% 5. Compressibility in accumulation term
%For component 1
CompAcc = x1(:,1).*dRhodP(:,1).*(pv/dt).*s(:,1) + x1(:,2).*dRhodP(:,2).*(pv/dt).*(1-s(:,1));
J1P = J1P + spdiags(-CompAcc,0,N,N);
%For component 2
CompAcc = (1-x1(:,1)).*dRhodP(:,1).*(pv/dt).*s(:,1) + (1-x1(:,2)).*dRhodP(:,2).*(pv/dt).*(1-s(:,1));
J2P = J2P + spdiags(-CompAcc,0,N,N);

%% Add capillarity
CapJwS = J2P * spdiags(dPc, 0, N, N);
J2S = J2S - CapJwS.*spdiags(1-x1(:,2),0,N,N);
J1S = J1S - CapJwS.*spdiags(x1(:,2),0,N,N);

%% Add wells
[J1P, J2P, J1S, J2S] = AddWellsToJacobian2(J1P, J2P, J1S, J2S, Inj, Prod, K, Mw, Mnw, dMw, dMnw, Status, Fluid);

%% Full Jacobian: put the 4 blocks together
J = [J1P, J1S; J2P, J2S];
end