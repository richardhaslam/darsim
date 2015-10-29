function [A, Tx, Ty] = AssemblePressureMatrix(Grid, K)
%Assemble Pressure Matrix

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
end