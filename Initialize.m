%Initialization the reservoir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Barnaby Fryer
%TU Delft
%Created: 3 June 2016 
%Last modified: 9 June 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Status] = Initialize(Fluid, Grid, K, FlashSettings)
%Initializes pressure, saturation, phase mole fraction, and total mole
%fraction in each cell of a 2-D reservoir based on either immiscible
%displacement, black oil, or compositional simulation.

switch (Fluid.Type)
    case('Immiscible')         %For the immiscible displacement case
        %Manually define pressure and saturation. x is simply [1 0] because no
        %interaction and z is simply equal to saturation
        
        Pinit = 0;                                  %HARD CODED initial pressure
        Swinit = .1;                                %HARD CODED initial saturation of phase 1
        
        Status.p = zeros(Grid.N, 1)*Pinit;          %Defines pressure vector
        Status.s = ones(Grid.N,2)*Swinit;           %Defines saturation vector for 2 phase
        Status.s(:,2) = (1-Swinit);                 %Defines second phase value
        
        Status.x1 = zeros(Grid.N, 2, 1);            %Defines phase mole fraction for 2 phase 2 component (each phase has its own column,
                                                    %values stored are for component 1)
        Status.x1(:,1) = Status.x1(:,1) + 1;        %Defines mole fraction of comp 1 in phase 1
        
        Status.z = Status.s;                        %For immiscible displacement total mole fraction is the saturation
        
    case('BlackOil')       %For the black oil model
        %Manually define pressure and saturation. x and z simply calculated
        
        Pinit = 1e5*linspace(50,10,Grid.N);                                  %HARD CODED initial pressure
        Swinit = .1;                                %HARD CODED initial saturation of phase 1
        
        Status.p = ones(Grid.N, 1).*Pinit';           %Defines pressure vector
        Status.s = ones(Grid.N,2)*Swinit;           %Defines saturation vector for 2 phase
        Status.s(:,2) = (1-Swinit);                 %Assigns second phase saturation
        Status.x1 = zeros(Grid.N,1);                %Predefines space for x1 vector
        Status.z = 'Initialize';                    %Sets z so that we do initialization in Innuer_Update.m
        
        %Do you need to send a dummy x!!!
        [Status] = Inner_Update(Status, Fluid, FlashSettings, Grid);   %Initializes reservoir using update function
        
        
    case('Compositional')  %For the compositional case
        %Manually define pressure and z. x and S are then calculated
        
        Pinit = 0;                                  %HARD CODED initial pressure
        zinit = .1;                                 %HARD CODED initial mole fraction of component 1
        
        Status.p = ones(Grid.N, 1)*Pinit;           %Defines pressure vector
        Status.s = 'Initialize';                    %Sets saturation so that we do initialization in Inner_Update.m
        Status.x1 = zeros(Grid.N,2);                %Predefines space for x1 vector
        Status.z = zeros(Grid.N, 2);                %Defines mole fraction vector for 2 component
        Status.z(:,1) = zinit;                      %Assigns first component mole fraction
        Status.z(:,2) = 1 - zinit;
        
        [Status] = Inner_Update(Status, Fluid, FlashSettings, Grid);   %Initializes reservoir using update function
        
end

%Removes column of second phase
Status.s = Status.s(:,1);
Status.z = Status.z(:,1);
Status.pc = ComputePc(Status.s, Fluid, K, Grid.por);
end

