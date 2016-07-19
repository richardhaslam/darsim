% Composition Update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 22 June 2016
%Last Modified: 22 June 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Composition Update
% Updates the composition and the saturation of the 2 phases.
% 2 options:
%   a. Immiscible: no compositional effects. Densisities can be different
%                  so mass fractions are updated.
%   b. Compositional fluid model: for Black oil or full compositional.

function [Status, delta, Converged] = CompositionUpdate(Status, Fluid, Grid, FlashSettings)              

switch(Fluid.Type)
    case('Immiscible')
        Status.s = min(Status.s,1);
        Status.s = max(Status.s,0);                                         %S does not get update in immiscible
        Status.z =  Compute_Z(Status.s,Status.x1, Status.rho);
        delta = 0;
        Converged = 1;
    otherwise    
        %% Flash loop
        [Status.z] =  Compute_Z(Status.s, Status.x1, Status.rho);          %Updates total mole fraction based on NR result
        
        Converged = 0;
        InnerCounter = 1;
        s_0 = Status.s;
        while Converged == 0 && InnerCounter < FlashSettings.MaxIt
            
            % 1. Stores old values
            s_old = Status.s;
            x_old = Status.x1;
            z_old = Status.z;
            
            % 2. Update Composition of the phases (Flash)
            [Status, SinglePhase] = Flash(Status, Fluid, Grid.Tres, FlashSettings.TolFlash);    
            
            %3. Update x and S based on components mass balance
            Status = ComputePhaseSaturation(Status, SinglePhase);
            
            %4. Compute new total mass fractions (z)
            [Status.z] =  Compute_Z(Status.s,Status.x1, Status.rho);
            
            %5.a Compute errors
            InnerError1 = norm(abs(Status.s(:,1) - s_old(:,1)), inf);   %Checks if this loop is converging
            InnerError2 = norm(abs(Status.x1(:,2) - x_old(:,2)), inf);
            InnerError3 = norm(abs(Status.z(:,1) - z_old(:,1)), inf);
            
            %5.b Check convergence
            if(InnerError1 < FlashSettings.TolInner && InnerError2 < FlashSettings.TolInner && InnerError3 < FlashSettings.TolInner)
                Converged = 1;
            end
            
            InnerCounter = InnerCounter + 1;
        end
        if Converged == 0               
            disp('Warning: the inner update has not converged! It may cause problems.');
        end
        delta = zeros(2*Grid.N, 1);
        delta(Grid.N+1:end) = Status.s - s_0;
end

end