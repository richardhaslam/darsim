function [Residual, T1, T2] = FIMResidualComp(Status0, Status, dt, Grid, K, Fluid, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod)

%Initialise local variables
p_old = Status0.p;
s_old = Status0.s;
x1_old = Status0.x1;
x2_old = 1-x1_old;
p = Status.p;
s = Status.s;
x1 = Status.x1;
x2 = 1 - x1;
pv = Grid.por*Grid.Volume;

%Density
[Rho, ~] = LinearDensity(p, Fluid.c, Fluid.rho);
[Rho_old, ~] = LinearDensity(p_old, Fluid.c, Fluid.rho);

%Accumulation term
A = speye(Grid.N)*pv/dt;
%Component 1
m1 = x1(:,1).*Rho(:,1).*s + x1(:,2).*Rho(:,2).*(1-s); 
m1_old = x1_old(:,1).*Rho_old(:,1).*s_old + x1_old(:,2).*Rho_old(:,2).*(1-s_old);
%Component 2
m2 = x2(:,1).*Rho(:,1).*s + x2(:,2).*Rho(:,2).*(1-s); 
m2_old = x2_old(:,1).*Rho_old(:,1).*s_old + x2_old(:,2).*Rho_old(:,2).*(1-s_old);

%Convective term
T1 = TransmissibilityMatrix (Grid, UpWindW, UpWindNw, Mw, Mnw, Rho, x1);
T2 = TransmissibilityMatrix (Grid, UpWindW, UpWindNw, Mw, Mnw, Rho, x2);

%Capillarity
%Matteo fixes this

%Gravity
G = ComputeGravityTerm(Grid.N);

%Source terms
[qw, qnw] = ComputeWells(Grid.N, Inj, Prod, K, Status, Mnw, Mw);

%% RESIDUAL
%Component 1
R1 = A * (m1 - m1_old)...
    + T1 * p...
    + G*s...
    - qw.*x1(:,1).*Rho(:,1) - qnw.*x1(:,2).*Rho(:,2);                      %Gravity and source terms
%Component 2
R2 = A * (m2 - m2_old)...
    + T2 * p...
    + G*s...
    - qw.*x2(:,1).*Rho(:,1) - qnw.*x2(:,2).*Rho(:,2);                      %Gravity and source terms

%Stick them together
Residual = [R1; R2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       EXTRA functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transmissibility
function T  = TransmissibilityMatrix(Grid, UpWindW, UpWindNw, Mw, Mnw, Rho, x)
%%%Transmissibility matrix construction
Tx = zeros(Grid.Nx+1, Grid.Ny);
Ty = zeros(Grid.Nx, Grid.Ny+1);

%Apply upwind operator
Mupx = UpWindW.x*(Mw .* Rho(:,1) .* x(:,1)) + UpWindNw.x*(Mnw .* Rho(:,2) .* x(:,2));
Mupy = UpWindW.y*(Mw .* Rho(:,1) .* x(:,1)) + UpWindNw.y*(Mnw .* Rho(:,2) .* x(:,2));
Mupx = reshape(Mupx, Grid.Nx, Grid.Ny);
Mupy = reshape(Mupy, Grid.Nx, Grid.Ny);
Tx(2:Grid.Nx,:)= Grid.Tx(2:Grid.Nx,:).*Mupx(1:Grid.Nx-1,:);
Ty(:,2:Grid.Ny)= Grid.Ty(:,2:Grid.Ny).*Mupy(:,1:Grid.Ny-1);

%Construct matrix 
x1 = reshape(Tx(1:Grid.Nx,:), Grid.N, 1);
x2 = reshape(Tx(2:Grid.Nx+1,:), Grid.N, 1);
y1 = reshape(Ty(:,1:Grid.Ny), Grid.N, 1);
y2 = reshape(Ty(:,2:Grid.Ny+1), Grid.N, 1);
DiagVecs = [-y2,-x2,y2+x2+y1+x1,-x1,-y1];
DiagIndx = [-Grid.Nx, -1, 0, 1, Grid.Nx];
T = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
end

%% Gravity
function G = ComputeGravityTerm(N)
G = speye(N)*0;
end

%% Wells
function [qw, qnw] = ComputeWells(N, Inj, Prod, K, Status, Mnw, Mw)
p = Status.p;
qw = zeros(N,1);
qnw = zeros(N,1);
%Injectors
for i=1:length(Inj)
    c = Inj(i).cells;
    switch (Inj(i).type)
        case('RateConstrained')
            qw(c) = Inj(i).q;
        case('PressureConstrained')
            qw(c) = Inj(i).Mw * Inj(i).PI .* K(c).* (Inj(i).p - p(c));
            qnw(c) = Inj(i).Mo * Inj(i).PI .* K(c) .* (Inj(i).p - p(c));
    end    
end
%Producers
for i=1:length(Prod)
    c = Prod(i).cells;
    switch (Prod(i).type)
        case('RateConstrained')
            qw(c) = Mw(c)./(Mw(c)+Mnw(c)).*Prod(i).q;
            qnw(c) = Prod(i).q - qw(c);
        case('PressureConstrained')
            qw(c) =   Mw(c).* Prod(i).PI .* K(c) .* (Prod(i).p - p(c));
            qnw(c) =  Mnw(c).* Prod(i).PI .* K(c).* (Prod(i).p - p(c));
    end
end
end