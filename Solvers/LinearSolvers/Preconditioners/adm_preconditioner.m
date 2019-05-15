% ADM Preconditioner 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_preconditioner < preconditioner
    properties
        R
        P
        Solver
    end
    methods
        function SetUp(obj, ProductionSystem, DiscretizationModel)
            % Choose where to coarsen and build ADM grid
            DiscretizationModel.SelectADMGrid(ProductionSystem);
            
            % Construct R & P based on ADM grid
            DiscretizationModel.BuildADMOperators();
            
            [obj.R, obj.P] = DiscretizationModel.OperatorsHandler.AssembleFullOperators(); 
        end
        function x = Solve(A, rhs)
             % Restrict system
            rhs_c = obj.R * rhs;
            A_c = obj.R * A * obj.P;
            
            % Solve Coarse System
            xc = obj.Solver.Solve(A_c, rhs_c);
            
            % Prolong to fine-scale resolution
            x = obj.P * xc;
        end
    end
end