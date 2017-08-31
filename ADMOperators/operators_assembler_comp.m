%  ADM Full operators assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 27 September 2016
%Last modified: 4 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_assembler_comp < operators_assembler    
    methods
        function [R, P] = Assemble(obj, ADMRest, ADMProl)
            % Build full operators
            R = blkdiag(ADMRest, ADMRest, ADMRest, ADMRest);
            P = blkdiag(ADMProl{1}, ADMProl{2}, ADMProl{2}, ADMProl{2});
        end
    end
end