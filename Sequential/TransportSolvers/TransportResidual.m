%Compute Transport Residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Residual, A, CapJac, df] = TransportResidual(snew, s0, q, pv, U, dt, Fluid, Grid, K)
%viscous fractional flow function
[fw, df] = ComputeFractionalFlow(snew, Fluid);

% viscous fluxes matrix
A = SaturationMatrix(Grid,U,q);      % Total fluxes matrix
% capillary fluxes
if isempty(Fluid.Pc)
    B = 0*speye(Grid.N);
    CapJac=zeros(Grid.N);
else
    B = CapFluxMatrix(Grid, reshape(snew, Grid.Nx, Grid.Ny), Fluid, K);
    CapJac = CapJacobian(Fluid, K, reshape(snew, Grid.Nx, Grid.Ny), Grid);
end
% Compute residual
Residual = pv/dt*(snew-s0)  - max(q,0) - A*fw - B*ones(Grid.N,1);
end

function B = CapFluxMatrix(Grid, S, Fluid, K)
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Grid.N;
%Builds a matrix that sums capillary fluxes
CapFlux.x = zeros(Nx+1, Ny);
CapFlux.y = zeros(Ny, Ny +1);
s = reshape(S, N, 1);
C = Fluid.Cances(s);
C = reshape(C, Nx, Ny);
CapFlux.x(2:Nx, :) = Grid.Ax*(K(1:Nx-1,:).*C(1:Nx-1, :) - K(2:Nx, :).*C(2:Nx, :))/Grid.dx;
CapFlux.y(:,2:Ny) = Grid.Ay*(K(:,1:Ny-1).*C(:,1:Ny-1) - K(:, 2:Ny).*C(:, 2:Ny))/Grid.dy;
%
x1 = reshape(CapFlux.x(1:Nx,:),N,1);  
y1 = reshape(CapFlux.y(:,1:Ny),N,1);
x2 = reshape(CapFlux.x(2:Nx+1,:),N,1);
y2 = reshape(CapFlux.y(:,2:Ny+1),N,1);

%assemble matrix
DiagVecs = -x1 - y1 + x2 + y2; % diagonal vectors
DiagIndx = 0; % diagonal index
B = spdiags(DiagVecs,DiagIndx,N,N);
end

function CapJac = CapJacobian(Fluid, K, s, Grid)
%constructs capillary pressure Jacobian matrix
Nx = Grid.Nx;
Ny = Grid.Ny;
N = Grid.N;
[Mw, Mo] = Mobilities(s, Fluid);
Mt = Mw + Mo;
[~, dPc] = ComputePc(s, Fluid);

f = Grid.Ax/Grid.dx*K.*Mw.*Mo./Mt.*dPc;
g = Grid.Ax/Grid.dx*K.*Mw.*Mo./Mt.*dPc;

CapJac = zeros(N);
 %inner cells
for i=2:Nx-1
    for j=2:Ny-1
        m = i + Nx*(j-1);
        CapJac(m, m-Nx) = g(i, j-1);
        CapJac(m, m+Nx) = g(i, j+1);
        CapJac(m, m-1) = g(i-1, j);
        CapJac(m, m+1) = g(i+1, j);
        CapJac(m,m) = -2*(f(i,j) + g(i,j));
    end
end
%left and right rows 
left = 1;
right = Nx;
for j=2:Ny-1
    r = right + Nx*(j-1);
    l = left + Nx*(j-1);
    %
    CapJac(l, l-Nx) = g(left, j-1);
    CapJac(l, l+Nx) = g(left, j+1);
    CapJac(l, l+1) = f(left + 1, j);
    CapJac(l,l) = -f(left,j) - 2*g(left,j);
    %
    CapJac(r, r-Nx) = g(right, j-1);
    CapJac(r, r+Nx) = g(right, j+1);
    CapJac(r, r+1) = f(right-1, j);
    CapJac(r,r) = -f(right, j) - 2*g(right, j);
end
%top and bottom rows
bottom = 1;
top = Ny;
for i=2:Nx-1
   t = i + Nx*(top-1);
   b = i + Nx*(bottom-1);
   %
   CapJac(b, b-1) = f(i-1, bottom);
   CapJac(b, b+1) = f(i+1, bottom);
   CapJac(b, b+Nx) = g(i, bottom+1);
   CapJac(b,b) = -2*f(i, bottom) - g(i, bottom);
   %
   CapJac(t, t-1) = f(i+1, top);
   CapJac(t, t+1) = f(i-1, top);
   CapJac(t, t-Nx) = g(i, top-1);
   CapJac(t, t) = -2*f(i, top) - g(i, top); 
end
%corners
% b-l
CapJac(1, 2) = f(2, 1);
CapJac(1, Nx) = g(1, 2);
CapJac(1, 1) = -f(1, 1) - g(1, 1);
% b-r
CapJac(Nx, Nx-1) = f(Nx-1, 1);
CapJac(Nx, 2*Nx) = g(Nx, 2);
CapJac(Nx, Nx) = -f(Nx, Nx) - g(Nx, Nx);
% t-l
CapJac(N-Nx+1, N-Nx+2) = f(2, Ny);
CapJac(N-Nx+1, N-2*Nx+1) = g(1, Ny-1);
CapJac(N-Nx+1, N-Nx+1) = -f(1,Ny) - g(1,Ny);
% t-r
CapJac(N, N-Nx) = f(Nx, Ny-1);
CapJac(N, N-1) = g(Nx-1, Ny);
CapJac(N, N) = -f(Nx, Ny) - g(Nx, Ny);
%Let's make it sparse
CapJac = sparse(CapJac);
end