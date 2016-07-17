%Initialization of the reservoir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Barnaby Fryer and Matteo Cusini
%TU Delft
%Created: 22 June 2016 
%Last modified: 22 June 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of the reservoir
% Initializes pressure, saturation and composition of the reservoir.

function [Status] = InitializeReservoir(Fluid, Grid, K, FlashSettings)
%P_init = ones(Grid.N, 1)*10e5;
P_init = linspace(1e5, 10e4, Grid.N);
z_init = ones(Grid.N,1)*0.0;

%Initialize objects
Status.s = zeros(Grid.N, 1);
Status.x1 = zeros(Grid.N, 2);
Status.z = zeros(Grid.N, 2);

SinglePhase.onlyvapor = zeros(Grid.N, 1);
SinglePhase.onlyliquid = zeros(Grid.N, 1);

%1. Assign initial pressure and total composition
Status.p = P_init';      
Status.z(:,1) = z_init;      
Status.z(:,2) = 1-z_init;

[Status.rho, ~] = LinearDensity(Status.p, Fluid.c, Fluid.rho);
switch (Fluid.Type)
    case('Immiscible')         
        %Composition in this case is fixed to be 1 and 0
        Status.x1(:,1) = ones(N, 1);           
        Status.x1(:,2) = zeros(N, 1);   
        
        SinglePhase.onlyvapor (Status.z(:,1) == 1) = 1;
        SinglePhase.onlyliquid (Status.z(:,2) == 1) = 1;
        
        Status = ComputePhaseSaturation(Status, SinglePhase);
    otherwise
       % 2. Update Composition of the phases (Flash)
       [Status, SinglePhase] = Flash(Status, Fluid, Grid.Tres, FlashSettings.TolFlash);    
            
       %3. Update x and S based on components mass balance
       Status = ComputePhaseSaturation(Status, SinglePhase);
end
Status.pc = ComputePc(Status.s, Fluid, reshape(K(1,:,:), Grid.N, 1), Grid.por);
end