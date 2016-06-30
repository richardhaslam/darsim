function [Residual, T1, T2, Tw] = FIMResidualComp(Status0, Status, dt, Grid, K, Mnw, Mw, UpWindNw, UpWindW, Inj, Prod)

%Initialise local variables
p_old = Status0.p;
s_old = Status0.s;
x1_old = Status0.x1;
x2_old = 1-x1_old;
p = Status.p;
s = Status.s;
x1 = Status.x1;
x2 = 1 - x1;
s2 = 1 -s;
s2_old = 1 -s_old;
Pc = Status.pc;
pv = Grid.por*Grid.Volume;
%Density
Rho = Status.rho; 
Rho_old = Status0.rho; 

%Accumulation term
A = speye(Grid.N)*pv/dt;
%Component 1
m1 = x1(:,1).*Rho(:,1).*s + x1(:,2).*Rho(:,2).*s2; 
m1_old = x1_old(:,1).*Rho_old(:,1).*s_old + x1_old(:,2).*Rho_old(:,2).*s2_old;
%Component 2
m2 = x2(:,1).*Rho(:,1).*s + x2(:,2).*Rho(:,2).*s2; 
m2_old = x2_old(:,1).*Rho_old(:,1).*s_old + x2_old(:,2).*Rho_old(:,2).*s2_old;

%Convective term
T1 = TransmissibilityMatrix (Grid, UpWindW, UpWindNw, Mw, Mnw, Rho, x1);
T2 = TransmissibilityMatrix (Grid, UpWindW, UpWindNw, Mw, Mnw, Rho, x2);

%Capillarity
Tw = TransmissibilityMatrix (Grid, UpWindW, UpWindNw, Mw, Mnw, Rho, [ones(Grid.N,1), zeros(Grid.N,1)]);

%Gravity
G = ComputeGravityTerm(Grid.N);

%Source terms
[q1, q2] = ComputeWells(Grid.N, Inj, Prod, K, Status, Mnw, Mw, Rho);

%% RESIDUAL
%Component 1
R1 = A * m1 - A * m1_old...    %Accumulation term
    + T1 * p...                %Convective term
    - x1(:,1).*(Tw*Pc)...           %Capillarity
    + G*s...                   %Gravity
    - q1;                      %Source terms
%Component 2
R2 = A * m2 - A * m2_old...
    + T2 * p...                %Convective term
    - x2(:,1).*(Tw*Pc)...           %Capillarity
    + G*s...                   %Gravity
    - q2;                      %Source terms

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
function [q1, q2] = ComputeWells(N, Inj, Prod, K, Status, Mnw, Mw, Rho)
p = Status.p;
q1 = zeros(N,1);
q2 = zeros(N,1);
x1 = Status.x1;
x2 = 1 - x1;
%Injectors
for i=1:length(Inj)
    c = Inj(i).cells;
    switch (Inj(i).type)
        case('RateConstrained')
            q1(c) = Inj(i).q * Inj(i).x1(1,1) * Inj(i).rho(1,1) + Inj(i).q * Inj(i).x1(1,2) * Inj(i).rho(1,2);
        case('PressureConstrained')
            q1(c) = Inj(i).Mw * Inj(i).PI .* K(c).* (Inj(i).p - p(c)) .* Inj(i).x1(1,1) * Inj(i).rho(1,1)...
                + Inj(i).Mo * Inj(i).PI .* K(c).* (Inj(i).p - p(c)) .* Inj(i).x1(1,2) * Inj(i).rho(1,2);
            q2(c) = Inj(i).Mw * Inj(i).PI .* K(c).* (Inj(i).p - p(c)) .* Inj(i).x2(1,1) * Inj(i).rho(1,1)...
                + Inj(i).Mo * Inj(i).PI .* K(c).* (Inj(i).p - p(c)) .* Inj(i).x2(1,2) * Inj(i).rho(1,2);
    end    
end
%Producers
for i=1:length(Prod)
    c = Prod(i).cells;
    switch (Prod(i).type)
        case('RateConstrained')
            q1(c) = (Mw(c)./(Mw(c)+Mnw(c)).*Prod(i).q)*x1(c,1)*Rho(c,1) + (Mnw(c)./(Mw(c)+Mnw(c)).*Prod(i).q)*x1(c,2)*Rho(c,2);
            q2(c) = (Mw(c)./(Mw(c)+Mnw(c)).*Prod(i).q)*x2(c,1)*Rho(c,1) + (Mnw(c)./(Mw(c)+Mnw(c)).*Prod(i).q)*x2(c,2)*Rho(c,2);
        case('PressureConstrained')
            q1(c) =   Mw(c).* Prod(i).PI .* K(c) .* (Prod(i).p - p(c)) .* x1(c,1) .* Rho(c,1)...
                + Mnw(c).* Prod(i).PI .* K(c) .* (Prod(i).p - p(c)) .* x1(c,2) .* Rho(c,2);
            q2(c) =  Mw(c).* Prod(i).PI .* K(c) .* (Prod(i).p - p(c)) .* x2(c,1) .* Rho(c,1)...
                + Mnw(c).* Prod(i).PI .* K(c) .* (Prod(i).p - p(c)) .* x2(c,2) .* Rho(c,2);
    end
end
end