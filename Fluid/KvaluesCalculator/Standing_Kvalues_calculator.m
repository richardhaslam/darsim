% Standing K-values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 September 2016
%Last modified: 15 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Standing_Kvalues_calculator < Kvalues_calculator
    properties (Constant)
        Const_T = 9/5; %From K to R
        Const_p = 0.000145037738; %From Pa to psi
    end
    properties
        F
        a
        c
        Num
    end
    methods
        function k = Compute(obj, p, T, Components)
            %% Compute K values
            % Convert units to use Standing's correlation for K values
            T = T*obj.Const_T;                    
            P = p*obj.Const_p;       
            Tb = zeros(1, 2);
            b = zeros(1, 2);
            obj.a = zeros(length(p), 1);
            obj.c = zeros(length(p), 1);
            obj.Num = zeros(length(p), 2);
            
            for i=1:2
                Tb(i) = Components(i).Tb*(9/5);
                b(i) = Components(i).b*(9/5);
            end
            
            obj.F = b .* ((1./Tb) - (1/T));                         %Finds F factor as per Standing 1979
            obj.a = 1.2 + 4.5 * 10^-4 * P + 15 * 10^-8 * P.^2;       %a coefficient
            obj.c = 0.89 - 1.7 * 10^-4 * P - 3.5 * 10^-8 * P.^2;     %c coefficient
            obj.Num(:,1) = 10.^(obj.a + obj.c * obj.F(1));
            obj.Num(:,2) = 10.^(obj.a + obj.c * obj.F(2));
            k(:, 1) = obj.Num(:,1)./P;                   %K1 as per Standing 1979
            k(:, 2) = obj.Num(:,2)./P;                   %K2 as per Standing 1979
        end
        function dk = DKvalDp(obj, p)
            % Convert units to use Standing's correlation for K values                    
            P = p*obj.Const_p;
            da =  4.5 * 10^-4 + 15 * 10^-8 * 2 * P;
            dc = - 1.7 * 10^-4 * P - 3.5 * 10^-8 * 2 * P;
            dexp(:,1) = da + obj.F(1) .* dc;
            dexp(:,2) = da + obj.F(2) .* dc;
            dNum = log(10)*dexp .* obj.Num; 
            dk(:,1) =  (P.*dNum(:, 1) - obj.Num(:, 1))./P.^2;
            dk(:,2) =  (P.*dNum(:, 2) - obj.Num(:, 2))./P.^2;
        end
    end
end