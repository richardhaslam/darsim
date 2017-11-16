%  ADM Full operators assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 27 September 2016
%Last modified: 7 November 2017
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
        function [R, P] = Assemble(obj, DiscretizationModel, ProductionSystem, Residual)
            % 1. Choose grid resolution
            if obj.VarIndex == 1
                DiscretizationModel.SelectADMGrid(ProductionSystem, Residual);
            end
            
            % 2. Build full operators
            R = DiscretizationModel.OperatorsHandler.ADMRest;
            P = DiscretizationModel.OperatorsHandler.ADMProl{obj.VarIndex};
            if obj.NumberOfEq > obj.VarIndex
                for i = obj.VarIndex+1:obj.NumberOfEq
                    R = blkdiag(R, DiscretizationModel.OperatorsHandler.ADMRest);
                    P = blkdiag(P, DiscretizationModel.OperatorsHandler.ADMProl{2});
                end
            end
        end
    end
end