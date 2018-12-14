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
                % The equation number 3 (Energy balance in solid phase, i.e., Rock Temperature) is only for reservoir
                if i==3 && any(strcmp(ProductionSystem.Reservoir.State.Properties.keys,'Tr'))
                    R = blkdiag(R, DiscretizationModel.OperatorsHandler.ADMRest_Reservoir);
                else
                    R = blkdiag(R, DiscretizationModel.OperatorsHandler.ADMRest);
                end
                P = blkdiag(P, DiscretizationModel.OperatorsHandler.ADMProl{i});
            end
        end
    end
end