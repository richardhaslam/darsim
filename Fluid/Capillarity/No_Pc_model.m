% No Pc model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 8 September 2016
%Last modified: 8 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef  No_Pc_model < capillary_pressure_model
    properties
    end
    methods 
        function Pc = ComputePc(obj, S)
            %Compute Pc and dPc analytically
            Pc = zeros(length(S), 1); 
        end
        function dPc = dPcdS(obj, S)
            dPc = zeros(length(S), 1); 
        end
    end
end
