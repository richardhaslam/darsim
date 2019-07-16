% ADM Full operators assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_assembler < handle
    properties
        NumberOfEq
    end
    methods
        function obj = operators_assembler(n_eq)
            obj.NumberOfEq = n_eq;
        end
    end
    methods (Abstract)
        obj = Assemble(obj);
    end
end

