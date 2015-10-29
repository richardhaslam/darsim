function [P, U, Wells, A, Ab, q]=PressureSolver(Grid, K, Inj, Prod)
%PRESSURE Solver

%1.Compute transmissibilities using harmonic average.
Nx = Grid.Nx; Ny = Grid.Ny; 
dx = Grid.dx; dy = Grid.dy;
N = Grid.N;
Ax = Grid.Ax;
Ay = Grid.Ay;

%Harmonic average of perm.
Kx=zeros(Nx+1,Ny);
Ky=zeros(Nx, Ny+1);
Kx(2:Nx,:)=2*K(1,1:Nx-1,:).*K(1,2:Nx,:)./(K(1,1:Nx-1,:)+K(1,2:Nx,:));
Ky(:,2:Ny)=2*K(2,:,1:Ny-1).*K(2,:,2:Ny)./(K(2,:,1:Ny-1)+K(2,:,2:Ny));

%Transmissibility
Tx=zeros(Nx+1, Ny);
Ty=zeros(Nx, Ny+1);
Tx(2:Nx,:)= Ax/dx.*Kx(2:Nx,:);
Ty(:,2:Ny)=Ay/dy.*Ky(:,2:Ny);

%Construct pressure matrix
x1=reshape(Tx(1:Nx,:),N,1);
x2=reshape(Tx(2:Nx+1,:),N,1);
y1=reshape(Ty(:,1:Ny),N,1);
y2=reshape(Ty(:,2:Ny+1),N,1);
DiagVecs = [-y2,-x2,y2+x2+y1+x1,-x1,-y1];
DiagIndx = [-Nx,-1,0,1,Nx];
A = spdiags(DiagVecs,DiagIndx,N,N);
Ab = A;
q=zeros(N,1);

%Modify A to consider 2 wells. Otherwise the matrix is
%singular.
a=Inj.x+(Inj.y-1)*Nx;
b=(Prod.y-1)*Nx+Prod.x;
A(a,a)=A(a,a)+Inj.PI*K(1,Inj.x,Inj.y);
A(b,b)=A(b,b)+Prod.PI*K(1,Prod.x,Prod.y);
q(a)=Inj.PI*K(1,Inj.x,Inj.y)*Inj.p;
q(b)=Prod.PI*K(1,Prod.x,Prod.y)*Prod.p;
%Solve for pressure
p = A\q;

%Compute total fluxes ([m^3/s])
P=reshape(p, Nx, Ny,1);
U.x = zeros(Nx+1,Ny,1);
U.y = zeros(Nx,Ny+1,1);
U.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)).*Tx(2:Nx,:);
U.y(:,2:Ny)  = (P(:,1:Ny-1)-P(:,2:Ny)).*Ty(:,2:Ny);

%Wells: fluxes [m^3/s]
Wells.Fluxes=zeros(Nx,Ny,1);
Wells.fw=zeros(Nx,Ny,1);
Wells.Fluxes(Inj.x,Inj.y)=Inj.PI*K(1,Inj.x,Inj.y)*(Inj.p-P(Inj.x,Inj.y));
Wells.Fluxes(Prod.x, Prod.y)=Prod.PI*K(1,Prod.x,Prod.y)*(Prod.p-P(Prod.x,Prod.y));
end