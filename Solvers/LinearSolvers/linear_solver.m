% Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef linear_solver < handle
    properties
        CondNumber
        %Preconditioner
    end
    methods
        function SetUp(obj, ProductionSystem, DiscretizationModel)
            %obj.Preconditioner.Setup(ProductionSystem, DiscretizationModel);
        end
        function x = Solve(obj, A, rhs)
            %x = obj.Preconditioner.Solve(A, rhs);
            x = A\rhs;
        end
    end
end