%  ADM Full operators assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 27 September 2016
%Last modified: 7 November 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_assembler_fim < operators_assembler
    properties
        
    end
    methods
        function obj = operators_assembler_fim(n_eq)
            obj@operators_assembler(n_eq);
        end
        function [R, P] = Assemble(obj, DiscretizationModel, ProductionSystem)
            % 1. Choose grid resolution
            DiscretizationModel.SelectADMGrid(ProductionSystem);
            
            % 2. Build full operators
            R = DiscretizationModel.OperatorsHandler.ADMRest;
            P = DiscretizationModel.OperatorsHandler.ADMProl{1};
            for i = 2:obj.NumberOfEq
                R = blkdiag(R, DiscretizationModel.OperatorsHandler.ADMRest);
                P = blkdiag(P, DiscretizationModel.OperatorsHandler.ADMProl{2});
            end
        end
    end
end