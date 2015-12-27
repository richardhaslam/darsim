%Build FIM Jacobian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = BuildJacobian_split(Grid, K, Trx, Try, P, Mw, Mo, dMw, dMo, U, dt, Inj, Prod, UpWind)
%Build FIM Jacobian
Nx=Grid.Nx; Ny=Grid.Ny; 
N=Nx*Ny;
pv=Grid.Volume*Grid.por;
%Wells' indexes in vectors
a=Inj.x+(Inj.y-1)*Nx;
b=(Prod.y-1)*Nx+Prod.x;

% BUILD FIM JACOBIAN BLOCK BY BLOCK

%1. Rpp Block
%Transmissibilitywith upwind oil mobility
Tx=zeros(Nx+1, Ny);
Ty=zeros(Nx, Ny+1);
Mupxo = UpWind.x*Mo;
Mupyo = UpWind.y*Mo;
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
Jpp = spdiags(DiagVecs,DiagIndx,N,N);
%Wells: producer only
Jpp(a,a)=Jop(a,a)+Inj.PI*(K(1,Inj.x,Inj.y)*K(2, Inj.x, Inj.y))^0.5*Inj.Mo;
Jpp(b,b)=Jop(b,b)+Prod.PI*(K(1,Prod.x,Prod.y)*K(2, Prod.x, Prod.y))^0.5*Mo(b);

%2. Rps Block
dMupxo = UpWind.x*dMo;
dMupyo = UpWind.y*dMo;
%Construct Jos block
x1=min(reshape(U.x(1:Nx,:),N,1),0).*dMupxo;
x2=max(reshape(U.x(2:Nx+1,:),N,1),0).*dMupxo;
y1=min(reshape(U.y(:,1:Ny),N,1),0).*dMupyo;
y2=max(reshape(U.y(:,2:Ny+1),N,1),0).*dMupyo;
v=ones(N,1)*pv/dt;
DiagVecs = [-y2,-x2,y2+x2-y1-x1-v,x1,y1];
DiagIndx = [-Nx,-1,0,1,Nx];
Jps = spdiags(DiagVecs,DiagIndx,N,N);
%Wells: Producer only
%Jps(a,a)=Jos(a,a)-Inj.PI*(K(1,Inj.x,Inj.y)*K(2, Inj.x, Inj.y))^0.5*(Inj.p-P(Inj.x,Inj.y))*Inj.dMo;
Jps(b,b)=Jos(b,b)-Prod.PI*(K(1,Prod.x,Prod.y)*K(2, Prod.x, Prod.y))^0.5*(Prod.p-P(Prod.x,Prod.y))*dMo(b);

%3. Rsp Block
%Transmissibility with upwind water mobility
Fx=zeros(Nx+1, Ny);
Fy=zeros(Nx, Ny+1);
Mupxw = UpWind.x*Mw;
Mupyw = UpWind.y*Mw;
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
Jsp = spdiags(DiagVecs,DiagIndx,N,N);
%Wells: Inj and prod
Jsp(a,a)=Jwp(a,a)+Inj.PI*(K(1,Inj.x,Inj.y)*K(2, Inj.x, Inj.y))^0.5*Inj.Mw;
Jsp(b,b)=Jwp(b,b)+Prod.PI*(K(1,Prod.x,Prod.y)*K(2, Prod.x, Prod.y))^0.5*Mw(b);

%4. Rss Block
dMupxw = UpWind.x*dMw;
dMupyw = UpWind.y*dMw;
%Construct Jws block
x1=min(reshape(U.x(1:Nx,:),N,1),0).*dMupxw;
x2=max(reshape(U.x(2:Nx+1,:),N,1),0).*dMupxw;
y1=min(reshape(U.y(:,1:Ny),N,1),0).*dMupyw;
y2=max(reshape(U.y(:,2:Ny+1),N,1),0).*dMupyw;
v=ones(N,1)*pv/dt;
DiagVecs = [-y2,-x2,y2+x2-y1-x1+v,x1,y1];
DiagIndx = [-Nx,-1,0,1,Nx];
Jss = spdiags(DiagVecs,DiagIndx,N,N);
%Wells: Inj and Prod
%Jss(a,a)=Jws(a,a)-Inj.PI*(K(1,Inj.x,Inj.y)*K(2, Inj.x, Inj.y))^0.5*(Inj.p-P(Inj.x,Inj.y))*Inj.dMw; 
Jss(b,b)=Jws(b,b)-Prod.PI*(K(1,Prod.x,Prod.y)*K(2, Prod.x, Prod.y))^0.5*(Prod.p-P(Prod.x,Prod.y))*dMw(b);

% Full Jacobian: put the 4 blocks together
J = [Jpp, Jps; Jsp, Jss];
end