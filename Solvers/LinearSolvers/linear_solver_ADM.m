% ADM Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef linear_solver_ADM < linear_solver_MMs
    properties
        Name
        Tol
        Maxit
        Iter
        R
        P
        Smooth = 0
    end
    methods
        function obj = linear_solver_ADM(name, tol, maxit)
           obj@linear_solver_MMs(name, tol, maxit)
        end
        function SetUp(obj, ProductionSystem, DiscretizationModel)
            % Choose where to coarsen and build ADM grid
            DiscretizationModel.SelectADMGrid(ProductionSystem);
            
            % Construct R & P based on ADM grid
            DiscretizationModel.BuildADMOperators();
            
            [obj.R, obj.P] = DiscretizationModel.AssembleFullOperators(); 
        end
    end
end