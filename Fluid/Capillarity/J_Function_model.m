% Capillary Pressure model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 8 September 2016
%Last modified: 8 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef  J_Function_model < capillary_pressure_model
    properties
        PorPermTerm
    end
    properties (Constant)
        alpha = 4.361e-2;
    end
    methods 
        function obj = J_Function_model(ProductionSystem)
            K = ProductionSystem.Reservoir.K(:,1);
            por = ProductionSystem.Reservoir.Por;
            obj.PorPermTerm = (por./K).^(0.5);
        end
        function Pc = ComputePc(obj, S)
            % Capillary Pressure
            
            %J-leverett curve
            
            J = 0.1.*(S).^(-0.5);
            
            %Compute Pc and dPc analytically
            Pc = obj.alpha .* obj.PorPermTerm.*J; 
        end
        function dPc = dPcdS(obj, S)
            % Derivative
            dJ = - 0.1*0.5*(S).^(-1.5);
            dPc = obj.alpha .* obj.PorPermTerm .* dJ;
        end
    end
end
