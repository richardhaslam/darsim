% Operators handler base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_handler < handle
    properties
        ProlongationBuilders
    end
    methods
        function obj = operators_handler(cf)
            obj.ProlongationBuilders = prolongation_builder.empty;
        end
        function AddProlongationBuilder(obj, prolongationbuilder, index)
            obj.ProlongationBuilders(index,1) = prolongationbuilder;
        end
    end
end