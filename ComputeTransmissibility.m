function [Tx, Ty] = ComputeTransmissibility(Grid, K)
%Transmissibilities
%%%%%%%%%%%%%%%%%%%%%%%
%Give a Grid Harmonic average of the permeability is computed and
%assigned to the interfaces. Then transmissibilities are computed.
%%%%%%%%%%%%%%%%%%%%%%%

Nx = Grid.Nx; 
Ny = Grid.Ny;
dx = Grid.dx;
dy = Grid.dy;
Ax = Grid.Ax;
Ay = Grid.Ay;

%Harmonic average of rock permeability.
Kx=zeros(Nx+1,Ny);
Ky=zeros(Nx, Ny+1);
Kx(2:Nx,:)=2*K(1,1:Nx-1,:).*K(1,2:Nx,:)./(K(1,1:Nx-1,:)+K(1,2:Nx,:));
Ky(:,2:Ny)=2*K(2,:,1:Ny-1).*K(2,:,2:Ny)./(K(2,:,1:Ny-1)+K(2,:,2:Ny));

%Transmissibility
Tx=zeros(Nx+1, Ny);
Ty=zeros(Nx, Ny+1);
Tx(2:Nx,:)= Ax/dx.*Kx(2:Nx,:);
Ty(:,2:Ny)=Ay/dy.*Ky(:,2:Ny);
end