%Flash calculator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 22 June 2016
%Last Modified: 23 June 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Flash calculator
% Updates the composition of the 2 phases. Two options:
%   1. Black oil: it computes an Rs that is assigned to x12 
%   2. Compositional:
%        This solves for the mole fractions in each phase for a two phase multi
%        component system. This uses the Standing (1979) correlation and is really
%        only valid for p < 1000 psia & T < 200 F.
%        Input: - Status.p [Pa]: pressure
%               - Status.z []: total mole/mass fraction of each component  
%               - Tres [K]: reservoir temperature
%               - Fluid.Comp.b [K]: the slope of the straight line connecting the critical point and the atmospheric boiling point on a log vapor pressure vs 1/T plot
%               - Fluid.Comp.Tb [K]: boiling point temperature of each component at the reference pressure (patm = 10^5 Pa)                        
        

function [Status, SinglePhase] = Flash(Status, Fluid, Tres, tol)

%1. Initialize local variables
p = Status.p;
z = Status.z;
N = length(p);
x = zeros(N, 2);

SinglePhase.onlyliquid = zeros(N, 1);
SinglePhase.onlyvapor = zeros(N, 1);

switch Fluid.Type
    case('BlackOil')
        %% - Dimensionless pressure!
        Pdim = p/Fluid.Pref;                  
        
        %% - Solve for x's
        x(:,2) = 1 - (800./(800 + 100*(0.2*Pdim(:,1) + 0.2)));  %More pressure less comp 1 (oil) in phase 1 (oil) (actually more gas pushed in really)
        x(:,1) = 1 - 0;                                         %light component is all in gas
        
        %Recognize single phase cells and fix their xs to be equal to z
        SinglePhase.onlyliquid(x(:, 2) >= z(:,1)) = 1;
        x(SinglePhase.onlyliquid == 1, 2) = z(SinglePhase.onlyliquid == 1, 1);
    case('Compositional')
        
        % 1. Convert units to use Standing's correlation for K values
        T = Tres*(9/5);                    %From K to R
        b = Fluid.Comp.b*(9/5);            %From K to R
        Tb = Fluid.Comp.Tb*(9/5);          %From K to R
        P = Status.p*0.000145037738;       %From Pa to psi 
        
        % 2. Compute K values
        F = b .* ((1./Tb) - (1/T));                         %Finds F factor as per Standing 1979
        a = 1.2 + 4.5 * 10^-4 * P + 15 * 10^-8 * P.^2;       %a coefficient
        c = 0.89 - 1.7 * 10^-4 * P - 3.5 * 10^-8 * P.^2;     %c coefficient
        k(:, 1) = (10.^(a + c * F(1)))./P;                   %K1 as per Standing 1979
        k(:, 2) = (10.^(a + c * F(2)))./P;                   %K2 as per Standing 1979
              
        
        % 3 Chek if we are in 2 phase region 
        
        % 3.a: checking if it 's all liquid: checks if mix is below bubble
        % point
        BubCheck = z .* k;                  
        BubCheck = sum(BubCheck, 2);
        
        x(BubCheck < 1, 2) = z (BubCheck < 1, 1);
        x(BubCheck < 1, 1) = 1;                               % This is to avoid having singular Jacobian matrix.
        SinglePhase.onlyliquid(BubCheck < 1) = 1;
        
        % 3.b: checking if it 's all vapor: checks if mix is above dew
        % point
        DewCheck = z ./ k;                  
        DewCheck = sum(DewCheck, 2);
        x(DewCheck < 1, 1) = z (DewCheck < 1, 1);
        x(DewCheck < 1, 2) = 0;
        SinglePhase.onlyvapor(DewCheck < 1) = 1;
                
        % 4. Actual Flash: solves for fv (vapor fraction)
        TwoPhase = ones(N, 1);
        TwoPhase(SinglePhase.onlyliquid == 1) = 0;
        TwoPhase(SinglePhase.onlyvapor == 1) = 0;
        
        alpha = 0.5; 
        
        %Initilaize variables
        fv = .5 * ones(N,1);   % 50-50 split as inital guess 
        %Single phase cells do not need to flash
        fv(SinglePhase.onlyvapor == 1) = 1;
        fv(SinglePhase.onlyliquid == 1)= 0;
        
        % Find fv with the tangent method
        converged = 0;
        while ~converged && alpha > 0.1
            itCounter = 0;
            while itCounter < 200 && ~converged
                %Finds hi for each component
                hi(:,1) = (z(:,1) .* k(:,1)) ./ (fv .* (k(:,1) - 1) + 1);             
                hi(:,2) = (z(:,2) .* k(:,2)) ./ (fv .* (k(:,2) - 1) + 1);            
                %Finds the derivative of hi for each component
                dhi(:, 1) = (z(:,1) .* (k(:,1) - 1).^2) ./ ((fv .* (k(:,1) - 1) + 1).^2);
                dhi(:, 2) = (z(:,2) .* (k(:,2) - 1).^2) ./ ((fv .* (k(:,2) - 1) + 1).^2);
                h = sum(hi, 2) - 1; 
                dh = - sum(dhi, 2); 
                
                %Update fv 
                h(TwoPhase == 0) = 0;
                fvnew = alpha * (-h ./ dh) + fv;
                                   
                fv = fvnew;
                if norm(h, inf) < tol
                    converged = 1;
                end
                itCounter = itCounter + 1;
            end
            alpha = alpha/2;
        end
        if ~converged 
            disp('Warning: The flash calculator could not converge.');
            return
        end
        
        %5. Solve for x's and y's
        x(TwoPhase == 1, 2) = z(TwoPhase == 1, 1) ./ (fv(TwoPhase == 1, 1) .* (k(TwoPhase == 1, 1) - 1) + 1);    %Solves for mole fractions in liquid phase
        x(TwoPhase == 1, 1) = k(TwoPhase == 1, 1) .* x(TwoPhase == 1, 2);                      %Solves for mole fractions in gas phase
end
Status.x1 = x;
end