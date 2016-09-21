% Gravity model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 20 September 2016
%Last modified: 20 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef gravity_model < handle
    properties (Constant)
        g = 9.806 % gravitational constant []
    end
    properties
        Depth % 
    end
    methods
        function ComputeBouyancyTerm()
            
        end
    end
end

