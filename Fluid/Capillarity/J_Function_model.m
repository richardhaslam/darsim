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
        WettingPhaseIndex
    end
    properties (Constant)
        alpha = 4.361e-2;
    end
    methods 
        function obj = J_function_model(wetting)
            obj.WettingPhaseIndex = wetting;
        end
        function Pc = ComputePc(obj, S, Phases)
            %J-leverett curve
            S = (S - Fluid.sr(2))./(1 - Fluid.sr(2)) + 0.1;
            J = 0.1.*(S).^(-0.5);
            dJ = - 0.1*0.5*(S).^(-1.5);
            %Define parameters
            
            %Compute Pc and dPc analytically
            Pc = obj.alpha .* (por./K).^(0.5).*J;
            dPc = obj.alpha .* (por./K).^(0.5).* dJ;  
        end
        function dPcdS(obj)
        end
    end
end
