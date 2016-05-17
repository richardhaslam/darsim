%Pressure Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, U, Pc, Wells, Inj, Prod, A, Ab, q]=PressureSolver(Grid, Inj, Prod, Fluid, S, K)
%PRESSURE Solver

%1.Compute transmissibilities using harmonic average.
Nx = Grid.Nx; Ny = Grid.Ny; 
N = Grid.N;

%Transmissibility
% Effective permeability
[Mw, Mo] = Mobilities(S, Fluid);
Mt = Mw+Mo;   %total mobility
Kt = zeros(2, Grid.Nx, Grid.Ny);
Kt(1,:,:) = reshape(Mt, 1, Grid.Nx, Grid.Ny).*K(1,:,:);		% x-direction
Kt(2,:,:) = reshape(Mt, 1, Grid.Nx, Grid.Ny).*K(2,:,:);		% y-direction
Ktvector = reshape(Kt(1,:,:), N, 1);
Kvector = reshape(K(1,:,:), N, 1);
[Tx, Ty] = ComputeTransmissibility(Grid, Kt);

%Construct pressure matrix
A = AssemblePressureMatrix(Tx, Ty, Nx, Ny);
Ab = A;
q = zeros(N,1);

%Add capillary term to the right-hand side
Pc = zeros(Nx, Ny);
Ucap.x = zeros(Nx+1,Ny);
Ucap.y = zeros(Nx, Ny+1);
if ~isempty(Fluid.Pc)
    Kw = zeros(2, Grid.Nx, Grid.Ny);
    Kw(1,:,:) = reshape(Mw, 1, Grid.Nx, Grid.Ny).*K(1,:,:);		% x-direction
    Kw(2,:,:) = reshape(Mw, 1, Grid.Nx, Grid.Ny).*K(2,:,:);		% y-direction
    [q, Pc, Ucap] = AddPcToPressureSystem(q, S, Fluid, Kw, K(1,:,:), Grid);
end

%% Add Wells
%Injectors
for i=1:length(Inj)
    a = Inj(i).cells;
    switch (Inj(i).type)
        case('PressureConstrained')
            for ii=1:length(a)
                A(a(ii),a(ii)) = A(a(ii),a(ii)) + Inj(i).PI*Kvector(a(ii)).*Inj(i).Mw;
                q(a(ii)) = q(a(ii)) + Inj(i).PI.*Kvector(a(ii)).*Inj(i).Mw.*Inj(i).p;
            end
        case('RateConstrained')
            for ii=1:length(a)
                q(a(ii)) = q(a(ii)) + Inj(i).q;
            end
    end
end
%Producers
for i=1:length(Prod)
    a = Prod(i).cells;
    switch(Prod(i).type)
        case('PressureConstrained')
            for ii=1:length(a)
                A(a(ii),a(ii)) = A(a(ii),a(ii)) + Prod(i).PI*Ktvector(a(ii));
                q(a(ii)) = q(a(ii)) + Prod(i).PI*Ktvector(a(ii)).*Prod(i).p;
            end
        case('RateConstrained')
            for ii=1:length(a)
                q(a(ii)) = q(a(ii)) + Prod(i).q;
            end
    end
end

%Solve for pressure
p = A\q;

%Compute total fluxes ([m^3/s])
P = reshape(p, Nx, Ny,1);
U.x = zeros(Nx+1,Ny,1);
U.y = zeros(Nx,Ny+1,1);
U.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)).*Tx(2:Nx,:) - Ucap.x(2:Nx,:);
U.y(:,2:Ny) = (P(:,1:Ny-1)-P(:,2:Ny)).*Ty(:,2:Ny) - Ucap.y(:,2:Ny);

%% Wells: fluxes [m^3/s]
Fluxes = zeros(N,1);
%Injectors
for i=1:length(Inj)
    a = Inj(i).cells;
    switch(Inj(i).type)
        case('PressureConstrained')
            Fluxes(a) = Fluxes(a) + Inj(i).PI.* Kvector(a).*Inj(i).Mw.*(Inj(i).p-p(a));
        case('RateConstrained')
            Fluxes(a) = Fluxes(a) + Inj(i).q;
            Inj(i).p = Inj(i).q./(Kvector(a).*Inj(i).Mw) + p(a);
    end
end
%Producers
for i=1:length(Prod)
    a = Prod(i).cells;
    switch(Prod(i).type)
        case('PressureConstrained')
            Fluxes(a) = Fluxes(a) + Prod(i).PI.* Ktvector(a).*(Prod(i).p-p(a));
        case('RateConstrained')
            Fluxes(a) = Fluxes(a) + Prod(i).q;
            Prod(i).p = Prod(i).q./(Ktvector(a)) + p(a);
    end
end
Wells.Fluxes = reshape(Fluxes, Nx, Ny);
end

