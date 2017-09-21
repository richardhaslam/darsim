% ADM Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 21 September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef linear_solver_ADM < linear_solver_MMs
    properties
        OperatorsAssembler
    end
    methods
        function obj = linear_solver_ADM(name, tol, maxit)
           obj@linear_solver_MMs(name, tol, maxit)
        end
        function SetUp(obj, DiscretizationModel)
            % Get ADM Operators
            [obj.R, obj.P] = obj.OperatorsAssembler.Assemble(DiscretizationModel.OperatorsHandler.ADMRest, DiscretizationModel.OperatorsHandler.ADMProl); 
        end
    end
end