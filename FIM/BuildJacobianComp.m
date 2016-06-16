%Build Jacobean matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 13 June 2016
%Last Modified: 13 June 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = BuildJacobianComp(Grid, K, TMatrix1, TMatrix2, Status, Mw, Mnw, dMw, dMnw, Rho, dRho, Uw, Unw, dPc, dt, Inj, Prod, UpWindW, UpWindNw)
%Build FIM Jacobian
Nx = Grid.Nx; 
Ny = Grid.Ny; 
N = Grid.N;
pv = Grid.Volume*Grid.por;
x1 = Status.x1;
x2 = 1 - x1;
s = Status.s;

% BUILD FIM JACOBIAN BLOCK BY BLOCK

%% 1 Component 1 pressure block

% 1.a: divergence
J1p = TMatrix1;

% 1.b: compressibility part
dMupxW = UpWindW.x*(Mw .* x1(:,1) .* dRho(:,1)); 
dMupyW = UpWindW.y*(Mw .* x1(:,1) .* dRho(:,1)); 
dMupxNw = UpWindNw.x*(Mnw .* x1(:,2) .* dRho(:,2));
dMupyNw = UpWindNw.y*(Mnw .* x1(:,2) .* dRho(:,2));

vecX1 = min(reshape(Uw.x(1:Nx,:),N,1),0).*dMupxW + min(reshape(Unw.x(1:Nx,:),N,1),0).*dMupxNw;
vecX2 = max(reshape(Uw.x(2:Nx+1,:),N,1),0).*dMupxW + max(reshape(Unw.x(2:Nx+1,:),N,1),0).*dMupxNw;
vecY1 = min(reshape(Uw.y(:,1:Ny),N,1),0).*dMupyW + min(reshape(Unw.y(:,1:Ny),N,1),0).*dMupyNw;
vecY2 = max(reshape(Uw.y(:,2:Ny+1),N,1),0).*dMupyW + max(reshape(Unw.y(:,2:Ny+1),N,1),0).*dMupyNw; 
acc = pv/dt .* ( x1(:,1) .* dRho(:,1) .* s + x1(:,2) .* dRho(:,2) .* (1-s));
DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
DiagIndx = [-Nx, -1, 0, 1, Nx];
J1p = J1p + spdiags(DiagVecs,DiagIndx,N,N);

%% 2. Component 2 pressure block

% 2.a divergence
J2p = TMatrix2;

% 2.b: compressibility part 
dMupxW = UpWindW.x*(Mw .* x2(:,1) .* dRho(:,1)); 
dMupyW = UpWindW.y*(Mw .* x2(:,1) .* dRho(:,1)); 
dMupxNw = UpWindNw.x*(Mnw .* x2(:,2) .* dRho(:,2));
dMupyNw = UpWindNw.y*(Mnw .* x2(:,2) .* dRho(:,2));

vecX1 = min(reshape(Uw.x(1:Nx,:),N,1),0).*dMupxW + min(reshape(Unw.x(1:Nx,:),N,1),0).*dMupxNw;
vecX2 = max(reshape(Uw.x(2:Nx+1,:),N,1),0).*dMupxW + max(reshape(Unw.x(2:Nx+1,:),N,1),0).*dMupxNw;
vecY1 = min(reshape(Uw.y(:,1:Ny),N,1),0).*dMupyW + min(reshape(Unw.y(:,1:Ny),N,1),0).*dMupyNw;
vecY2 = max(reshape(Uw.y(:,2:Ny+1),N,1),0).*dMupyW + max(reshape(Unw.y(:,2:Ny+1),N,1),0).*dMupyNw; 
acc = pv/dt .* (x2(:,1) .* dRho(:,1) .* s + x2(:,2) .* dRho(:,2) .* (1-s));
DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
DiagIndx = [-Nx, -1, 0, 1, Nx];
J2p = J2p + spdiags(DiagVecs,DiagIndx,N,N);

%% 3. Component 1 saturation block
dMupxW = UpWindW.x*(dMw .* x1(:,1) .* Rho(:,1)); 
dMupyW = UpWindW.y*(dMw .* x1(:,1) .* Rho(:,1)); 
dMupxNw = UpWindNw.x*(dMnw .* x1(:,2) .* Rho(:,2));
dMupyNw = UpWindNw.y*(dMnw .* x1(:,2) .* Rho(:,2));

vecX1 = min(reshape(Uw.x(1:Nx,:),N,1),0).*dMupxW + min(reshape(Unw.x(1:Nx,:),N,1),0).*dMupxNw;
vecX2 = max(reshape(Uw.x(2:Nx+1,:),N,1),0).*dMupxW + max(reshape(Unw.x(2:Nx+1,:),N,1),0).*dMupxNw;
vecY1 = min(reshape(Uw.y(:,1:Ny),N,1),0).*dMupyW + min(reshape(Unw.y(:,1:Ny),N,1),0).*dMupyNw;
vecY2 = max(reshape(Uw.y(:,2:Ny+1),N,1),0).*dMupyW + max(reshape(Unw.y(:,2:Ny+1),N,1),0).*dMupyNw; 
acc = pv/dt .* (x1(:,1) .* Rho(:,1) - x1(:,2) .* Rho(:,2));
DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
DiagIndx = [-Nx, -1, 0, 1, Nx];
J1S = spdiags(DiagVecs,DiagIndx,N,N);

%% 4. Component 2 saturation block
dMupxW = UpWindW.x*(dMw .* x2(:,1) .* Rho(:,1)); 
dMupyW = UpWindW.y*(dMw .* x2(:,1) .* Rho(:,1)); 
dMupxNw = UpWindNw.x*(dMnw .* x2(:,2) .* Rho(:,2));
dMupyNw = UpWindNw.y*(dMnw .* x2(:,2) .* Rho(:,2));

vecX1 = min(reshape(Uw.x(1:Nx,:),N,1),0).*dMupxW + min(reshape(Unw.x(1:Nx,:),N,1),0).*dMupxNw;
vecX2 = max(reshape(Uw.x(2:Nx+1,:),N,1),0).*dMupxW + max(reshape(Unw.x(2:Nx+1,:),N,1),0).*dMupxNw;
vecY1 = min(reshape(Uw.y(:,1:Ny),N,1),0).*dMupyW + min(reshape(Unw.y(:,1:Ny),N,1),0).*dMupyNw;
vecY2 = max(reshape(Uw.y(:,2:Ny+1),N,1),0).*dMupyW + max(reshape(Unw.y(:,2:Ny+1),N,1),0).*dMupyNw; 
acc = pv/dt .* (x2(:,1) .* Rho(:,1) - x2(:,2) .* Rho(:,2));
DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
DiagIndx = [-Nx, -1, 0, 1, Nx];
J2S = spdiags(DiagVecs,DiagIndx,N,N);

%% 5.Add capillarity


%% Add wells to each block
[J1p, J2p, J1S, J2S] = AddWellsToJacobeanComp(J1p, J2p, J1S, J2S, Inj, Prod, K, Status, Rho, dRho, Mw, Mnw, dMw, dMnw);
%% Full Jacobean: combine the 4 blocks
J = [J1p, J1S; J2p, J2S];

end