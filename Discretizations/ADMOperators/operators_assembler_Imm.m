%  ADM Full operators assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 27 September 2016
%Last modified: 27 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_assembler_Imm < operators_assembler
    methods
        function [R, P] = Assemble(obj, ADMRest, ADMProlp, ADMProls)
            % Build full operators
            R = blkdiag(ADMRest, ADMRest);
            P = blkdiag(ADMProlp, ADMProls);
        end
    end
end