%Multilevel Multiscale Finite Volume method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, U] = MMsFVPressureSolver(Grid, Inj, Prod, K, Fluid, S, CoarseGrid, maxLevel)
%MULTILEVEL MsFV method
%% AsemblePressureMatrix
%Transmissibility
% Effective permeability
[Mw, Mo] = Mobilities(S, Fluid);
Mt = Mw+Mo;   %total mobility
Kt = zeros(2, Grid.Nx, Grid.Ny);
Kt(1,:,:) = reshape(Mt, 1, Grid.Nx, Grid.Ny).*K(1,:,:);		% x-direction
Kt(2,:,:) = reshape(Mt, 1, Grid.Nx, Grid.Ny).*K(2,:,:);		% y-direction
Kvector = reshape(Kt(1,:,:), Grid.N, 1);
[Tx, Ty] = ComputeTransmissibility(Grid, Kt);
[Ap] = AssemblePressureMatrix(Tx, Ty, Grid.Nx, Grid.Ny);
q = zeros(Grid.N,1);

%Add Wells
pc = zeros(Grid.N, 1);
for i=1:length(Inj)
    [Ap, q] = AddWell(Ap, q, Inj(i), Kvector, pc, Kvector);
end
for i=1:length(Prod)
    [Ap, q] = AddWell(Ap, q, Prod(i), Kvector, pc, Kvector);
end

%% Construct R and P
[CoarseGrid(1).MsR, CoarseGrid(1).MsP, CoarseGrid(1).C] = MSFVOperators(Grid, CoarseGrid(1), Ap, 1);
CoarseGrid(1).A_c = CoarseGrid(1).MsP'*Ap*CoarseGrid(1).MsP;
CoarseGrid(1).q_c = CoarseGrid(1).MsP'*(q - Ap*CoarseGrid(1).C*q);
for x = 2:maxLevel
    [CoarseGrid(x).MsR, CoarseGrid(x).MsP, CoarseGrid(x).C] = MSFVOperators(CoarseGrid(x-1), CoarseGrid(x), CoarseGrid(x-1).A_c, x);
    CoarseGrid(x).A_c = CoarseGrid(x).MsP'*CoarseGrid(x-1).A_c*CoarseGrid(x).MsP;
    CoarseGrid(x).q_c = CoarseGrid(x).MsP'*(CoarseGrid(x-1).q_c - CoarseGrid(x-1).A_c*CoarseGrid(x).C*CoarseGrid(x-1).q_c);
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
P = reshape(p_f, Grid.Nx, Grid.Ny);

%% Reconstruct Velocity Field
%No reconstruction at the moment

%% Compute Fluxes
U.x = zeros(Grid.Nx+1,Grid.Ny,1);
U.y = zeros(Grid.Nx,Grid.Ny+1,1);
U.x(2:Grid.Nx,:) = (P(1:Grid.Nx-1,:)-P(2:Grid.Nx,:)).*Tx(2:Grid.Nx,:);
U.y(:,2:Grid.Ny)  = (P(:,1:Grid.Ny-1)-P(:,2:Grid.Ny)).*Ty(:,2:Grid.Ny);
end