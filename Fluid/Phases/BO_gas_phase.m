% BO gas phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 28 September 2016
%Last modified: 28 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_gas_phase < phase
    properties
        Bg
    end
    methods
        function rho = ComputeDensity(obj, Status, Components)
            rho = 1;
        end
        function drho = DrhoDp(obj, Status, Components)
            drho = 1;
        end
    end    
end