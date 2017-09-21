%  ADM Full operators assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 27 September 2016
%Last modified: 21 September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_assembler_fim < operators_assembler
    properties
        
    end
    methods
        function [R, P] = Assemble(obj, ADMRest, ADMProl)
            % Build full operators
            R = ADMRest;
            P = ADMProl{1};
            for i = 1:obj.NumberOfEq
                R = blkdiag(R, ADMRest);
                P = blkdiag(P, ADMProl{2});
            end
        end
    end
end