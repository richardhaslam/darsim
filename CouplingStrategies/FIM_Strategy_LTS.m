% FIM LTS coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FIM_Strategy_LTS < FIM_Strategy
    properties
        
    end
    
    methods
        function obj = untitled(inputArg1,inputArg2)
            obj.Property1 = inputArg1 + inputArg2;
        end
        function outputArg = method1(obj,inputArg)
           
            outputArg = obj.Property1 + inputArg;
        end
    end
end

