% ADM Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef linear_solver_ADM < linear_solver
    properties
        R
        P
    end
    methods
        function SetUp(obj, ProductionSystem, DiscretizationModel)
            % Choose where to coarsen and build ADM grid
            DiscretizationModel.SelectADMGrid(ProductionSystem);
            
            % Construct R & P based on ADM grid
            DiscretizationModel.BuildADMOperators();
            
            % Build full operators
            obj.R = blkdiag(DiscretizationModel.OperatorsHandler.ADMRest, DiscretizationModel.OperatorsHandler.ADMRest);
            obj.P = blkdiag(DiscretizationModel.OperatorsHandler.ADMProlp, DiscretizationModel.OperatorsHandler.ADMProls);
        end
        function xf = Solve(obj, A, rhs)
            % Restrict system
            rhs_c = obj.R * rhs;
            A_c = obj.R * A * obj.P;
            
            % Solve Coarse System
            x = A_c\rhs_c;
            
            % Prolong to fine-scale resolution
            xf = obj.P * x;
        end
    end
end