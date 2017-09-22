%  ADM Full operators assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 27 September 2016
%Last modified: 21 September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_assembler_seq < operators_assembler
    properties
        VarIndex
    end
    methods
        function obj = operators_assembler_seq(var_index, n_eq)
            obj@operators_assembler(n_eq);
            obj.VarIndex = var_index;
        end
        function [R, P] = Assemble(obj, ADMRest, ADMProl)
            % Build full operators
            R = ADMRest;
            P = ADMProl{obj.VarIndex};
            if obj.NumberOfEq > obj.VarIndex
                for i = obj.VarIndex+1:obj.NumberOfEq
                    R = blkdiag(R, ADMRest);
                    P = blkdiag(P, ADMProl{2});
                end
            end
        end
    end
end