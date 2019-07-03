classdef LTS_operators_assembler_seq < operators_assembler_seq
    methods
        function obj = LTS_operators_assembler_seq(var_index, n_eq)
            obj@operators_assembler_seq(var_index, n_eq);
        end
        function [R, P] = Assemble(obj, DiscretizationModel, ProductionSystem, Residual)
            % 1. Choose grid resolution
            if obj.VarIndex == 2
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
        
        function [R, P] = AssembleCoarse(obj, DiscretizationModel, ProductionSystem, Residual)
            % 1. Choose grid resolutio
            if obj.VarIndex == 2
                DiscretizationModel.SelectADMGridCoarse(ProductionSystem, Residual);
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