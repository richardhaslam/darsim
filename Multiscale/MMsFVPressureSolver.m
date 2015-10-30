function [P, U] = MMsFVPressureSolver(FineGrid, Inj, Prod, Kt, CoarseGrid, maxLevel)
%MULTILEVEL MsFV method
%% AsemblePressureMatrix
[Tx, Ty] = ComputeTransmissibility(FineGrid, Kt);
[Ap] = AssemblePressureMatrix(Tx, Ty);
q = zeros(FineGrid.N,1);

% Add Wells
a=Inj.x+(Inj.y-1)*FineGrid.Nx;
b=(Prod.y-1)*FineGrid.Nx+Prod.x;
Ap(a,a)=Ap(a,a)+Inj.PI*Kt(1,Inj.x,Inj.y);
Ap(b,b)=Ap(b,b)+Prod.PI*Kt(1,Prod.x,Prod.y);
q(a)=Inj.PI*Kt(1,Inj.x,Inj.y)*Inj.p;
q(b)=Prod.PI*Kt(1,Prod.x,Prod.y)*Prod.p;

%% Construct R and P
[CoarseGrid(1).MsR, CoarseGrid(1).MsP, CoarseGrid(1).C] = MSFVOperators(FineGrid, CoarseGrid(1), Ap, 1);
CoarseGrid(1).A_c = CoarseGrid(1).MsR*Ap*CoarseGrid(1).MsP;
CoarseGrid(1).q_c = CoarseGrid(1).MsR*(q - Ap*CoarseGrid(1).C*q);
for x = 2:maxLevel
    [CoarseGrid(x).MsR, CoarseGrid(x).MsP, CoarseGrid(x).C] = MSFVOperators(CoarseGrid(x-1), CoarseGrid(x), CoarseGrid(x-1).A_c, x);
    CoarseGrid(x).A_c = CoarseGrid(x).MsR*CoarseGrid(x-1).A_c*CoarseGrid(x).MsP;
    CoarseGrid(x).q_c = CoarseGrid(x).MsR*(CoarseGrid(x-1).q_c - CoarseGrid(x-1).A_c*CoarseGrid(x).C*CoarseGrid(x-1).q_c);
end
%% Solve
p_c = CoarseGrid(maxLevel).A_c\CoarseGrid(maxLevel).q_c;
%p_c = bicg(CoarseGrid(maxLevel).A_c, CoarseGrid(maxLevel).q_c, 1e-6,100,[],[],p_c);
for x = maxLevel:-1:2
    p_c = CoarseGrid(x).MsP*p_c + CoarseGrid(x).C*CoarseGrid(x-1).q_c;
    %p_c = bicg(CoarseGrid(x-1).A_c, CoarseGrid(x-1).q_c, 1e-6,100,[],[],p_c);
end
p_f = CoarseGrid(1).MsP*p_c + CoarseGrid(1).C*q;
%p_f = bicg(Ap, q, 1e-2, 100, [], [], p_f);
P = reshape(p_f, FineGrid.Nx, FineGrid.Ny);

%% Reconstruct Velocity Field
%No reconstruction at the moment

%% Compute Fluxes
U.x = zeros(FineGrid.Nx+1,FineGrid.Ny,1);
U.y = zeros(FineGrid.Nx,FineGrid.Ny+1,1);
U.x(2:FineGrid.Nx,:) = (P(1:FineGrid.Nx-1,:)-P(2:FineGrid.Nx,:)).*Tx(2:FineGrid.Nx,:);
U.y(:,2:FineGrid.Ny)  = (P(:,1:FineGrid.Ny-1)-P(:,2:FineGrid.Ny)).*Ty(:,2:FineGrid.Ny);
end