function [Status] = Inner_Update(Status, Fluid, FlashSettings, Grid)

InnerError = 1;                                             %Sets up error so it enters update loop
snew = zeros(Grid.N,2);                                     %Predefines memory for dummy saturation
[rho,~] = LinearDensity(Status.p, Fluid.c, Fluid.rho);        %Finds density in each cell for each phase

%Black Oil inner update
switch(Fluid.Type)
    case('BlackOil')
        InnerCounter = 0;
        Status.s(:,2) = 1 - Status.s(:,1);                      %Initializes counter
        if (strcmp(Status.z,'Initialize')==1)                   %Checks if we are doing initialization or have just come from NR loop
        else
            [Status.z] = Update_z(Status.s,Status.x1,rho);          %Updates total mole fraction based on NR result
        end
        
        while InnerError > FlashSettings.TolInner
            
            [Status.x1] = BO_Flash(Status.p);                   %Updates x with Black Oil Flash
            
            if (strcmp(Status.z,'Initialize')==1)               %Checks if we are doing initialization or have just come from NR loop
                snew = Status.s;                                %Updates dummy s that we can find z later and have zero error in error check
            else                                                %This is NOT initialization so update s and x
                snew(Status.x1(:,1) < Status.z(:,1),1) = 1;     %One phase S limit
                snew(Status.x1(:,1) < Status.z(:,1),2) = 0;
                [snew(Status.x1(:,1) > Status.z(:,1),:)] = Update_S(Status.x1(Status.x1(:,1) > Status.z(:,1),:),...
                    Status.z(Status.x1(:,1) > Status.z(:,1),:),rho(Status.x1(:,1) > Status.z(:,1),:));      %When two-phases
                
                Status.x1(Status.x1(:,1) < Status.z(:,1),1) = Status.z(Status.x1(:,1) < Status.z(:,1),1);           %One phase x limit
                Status.x1(Status.x1(:,2) > Status.z(:,2),2) = Status.z(Status.x1(:,2) > Status.z(:,2),2);
            end
            
            [Status.z] = Update_z(snew,Status.x1,rho);          %Updates z with new parameters. Different from my version
            
            InnerError = sum(abs(Status.s(:,1) - snew(:,1)));   %Checks if this loop is converging
            Status.s = snew;                                    %Makes dummy S the new real S
            InnerCounter = InnerCounter + 1;
        end
        
        %Compositional inner update
    case('Compositional')
        InnerCounter = 0;                                       %Initializes counter
        
        if (strcmp(Status.s,'Initialize')==1)                   %Checks if we are doing initialization or have just come from NR loop
            Status.x1 = zeros(Grid.N,2,1);                      %Predefines memory for x
        else
            Status.s(:,2) = 1 - Status.s(:,1);
            [Status.z] = Update_z(Status.s,Status.x1,rho);       %Updates total mole fraction based on NR result
        end
        
        while InnerError > FlashSettings.TolInner               %Checks error in loop
            for ii = 1:Grid.Nx*Grid.Ny
                [PhaseOneSplit,PhaseTwoSplit] = Comp_Calc(Status.p(ii,1),Grid.Tres,Fluid.Comp.b,Fluid.Comp.Tb,Status.z(ii,:),FlashSettings.TolFlash);    %Finds new phase mole fractions (x)
                Status.x1(ii,1) = PhaseOneSplit(1,1);           %Assigns phase split of phase 1
                Status.x1(ii,2) = PhaseTwoSplit(1,1);           %Assigns phase split of phase 2
            end
            [snew] = Update_S(Status.x1,Status.z,rho);          %Update s based on new x's
            [Status.z] = Update_z(snew,Status.x1,rho);          %Update z so that everything agrees
            InnerError = sum(abs(Status.s(:,1) - snew(:,1)));   %Find error of loop
            if (strcmp(Status.s,'Initialize')==1) 
                InnerError = 0;
            else
                Status.s = snew;                                    %Dummy s becomes s for next loop
            end
            InnerCounter = InnerCounter + 1;                    %Update counter
            if InnerCounter > FlashSettings.MaxIt               %If this loop iterates too much we stop it
                InnerError = 0;                                 %Allows exit of loop
            end
        end
        
        %Immiscible inner update
    case('immiscible')
        Status.s = Status.s;                                     %S does not get update in immiscible
        Status.x1 = zeros(Grid.N, 2, 1);                         %Phase split is always the same (every component in its own phase)
        Status.x1(:,1) = Status.x1(:,1) + 1;
        Status.z = Status.s;                                     %Total mole fraction is simply saturation in immiscible case
end
if (strcmp(Status.s,'Initialize')==1) 
    [Status.s] = Update_S(Status.x1,Status.z,rho);
else
    Status.s = Status.s(:,1);
end

end

