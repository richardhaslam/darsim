%  ADM Full operators assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_assembler_fim < operators_assembler
    properties
        
    end
    methods
        function obj = operators_assembler_fim(n_eq)
            obj@operators_assembler(n_eq);
        end
        function [R, P] = Assemble(obj, DiscretizationModel, ProductionSystem, Residual)
            % 1. Choose grid resolution and create operators
            DiscretizationModel.SelectADMGrid(ProductionSystem, Residual);
            
            % 2. Build full operators (the block diagonal matrices)
            R = DiscretizationModel.OperatorsHandler.ADMRest;
            P = DiscretizationModel.OperatorsHandler.ADMProl{1};
            for i = 2:obj.NumberOfEq
                R = blkdiag(R, DiscretizationModel.OperatorsHandler.ADMRest);
                P = blkdiag(P, DiscretizationModel.OperatorsHandler.ADMProl{2});
            end
        end
    end
end