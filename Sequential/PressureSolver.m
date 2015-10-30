function [P, U, Wells, A, Ab, q]=PressureSolver(Grid, K, Inj, Prod)
%PRESSURE Solver

%1.Compute transmissibilities using harmonic average.
Nx = Grid.Nx; Ny = Grid.Ny; 
N = Grid.N;

%Transmissibility
[Tx, Ty] = ComputeTransmissibility(Grid, K);

%Construct pressure matrix
A = AssemblePressureMatrix(Tx, Ty, Nx, Ny);
Ab = A;
q=zeros(N,1);

%Add Wells
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