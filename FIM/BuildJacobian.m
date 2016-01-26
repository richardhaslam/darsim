%Build FIM Jacobian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = BuildJacobian(Grid, K, Trx, Try, p, Mw, Mo, dMw, dMo, Uo, Uw, dt, Inj, Prod, UpWindO, UpWindW)
%Build FIM Jacobian
Nx=Grid.Nx; Ny=Grid.Ny; 
N=Nx*Ny;
pv=Grid.Volume*Grid.por;
Kvector = reshape(K(1, :, :), Grid.N, 1);

% BUILD FIM JACOBIAN BLOCK BY BLOCK

%1. Ro Pressure Block
%Transmissibilitywith upwind oil mobility
Tx=zeros(Nx+1, Ny);
Ty=zeros(Nx, Ny+1);
Mupxo = UpWindO.x*Mo;
Mupyo = UpWindO.y*Mo;
Mupxo = reshape(Mupxo, Nx, Ny);
Mupyo = reshape(Mupyo, Nx, Ny);
Tx(2:Nx,:)= Trx(2:Nx,:).*Mupxo(1:Nx-1,:);
Ty(:,2:Ny)= Try(:,2:Ny).*Mupyo(:,1:Ny-1);
%Construct Jop block
x1=reshape(Tx(1:Nx,:),N,1);
x2=reshape(Tx(2:Nx+1,:),N,1);
y1=reshape(Ty(:,1:Ny),N,1);
y2=reshape(Ty(:,2:Ny+1),N,1);
DiagVecs = [-y2,-x2,y2+x2+y1+x1,-x1,-y1];
DiagIndx = [-Nx,-1,0,1,Nx];
Jop = spdiags(DiagVecs,DiagIndx,N,N);

%2. Rw Pressure Block
%Transmissibility with upwind water mobility
Fx=zeros(Nx+1, Ny);
Fy=zeros(Nx, Ny+1);
Mupxw = UpWindW.x*Mw;
Mupyw = UpWindW.y*Mw;
Mupxw = reshape(Mupxw, Nx, Ny);
Mupyw = reshape(Mupyw, Nx, Ny);
Fx(2:Nx,:)= Trx(2:Nx,:).*Mupxw(1:Nx-1,:);
Fy(:,2:Ny)= Try(:,2:Ny).*Mupyw(:,1:Ny-1);
%Construct Jwp block
x1=reshape(Fx(1:Nx,:),N,1);
x2=reshape(Fx(2:Nx+1,:),N,1);
y1=reshape(Fy(:,1:Ny),N,1);
y2=reshape(Fy(:,2:Ny+1),N,1);
DiagVecs = [-y2,-x2,y2+x2+y1+x1,-x1,-y1];
DiagIndx = [-Nx,-1,0,1,Nx];
Jwp = spdiags(DiagVecs,DiagIndx,N,N);

%3. Ro Saturation Block
dMupxo = UpWindO.x*dMo;
dMupyo = UpWindO.y*dMo;
%Construct Jos block
x1=min(reshape(Uo.x(1:Nx,:),N,1),0).*dMupxo;
x2=max(reshape(Uo.x(2:Nx+1,:),N,1),0).*dMupxo;
y1=min(reshape(Uo.y(:,1:Ny),N,1),0).*dMupyo;
y2=max(reshape(Uo.y(:,2:Ny+1),N,1),0).*dMupyo;
v=ones(N,1)*pv/dt;
DiagVecs = [-y2,-x2,y2+x2-y1-x1-v,x1,y1];
DiagIndx = [-Nx,-1,0,1,Nx];
Jos = spdiags(DiagVecs,DiagIndx,N,N);

%4. Rw Saturation Block
dMupxw = UpWindW.x*dMw;
dMupyw = UpWindW.y*dMw;
%Construct Jws block
x1=min(reshape(Uw.x(1:Nx,:),N,1),0).*dMupxw;
x2=max(reshape(Uw.x(2:Nx+1,:),N,1),0).*dMupxw;
y1=min(reshape(Uw.y(:,1:Ny),N,1),0).*dMupyw;
y2=max(reshape(Uw.y(:,2:Ny+1),N,1),0).*dMupyw;
v=ones(N,1)*pv/dt;
DiagVecs = [-y2,-x2,y2+x2-y1-x1+v,x1,y1];
DiagIndx = [-Nx,-1,0,1,Nx];
Jws = spdiags(DiagVecs,DiagIndx,N,N);

% Add wells
[Jop, Jwp, Jos, Jws] = AddWellsToJacobian(Jop, Jwp, Jos, Jws, Inj, Prod, Kvector, p, Mw, Mo, dMw, dMo);

% Full Jacobian: put the 4 blocks together
J = [Jop, Jos; Jwp, Jws];
end