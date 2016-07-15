% Simulation class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 13 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Reservoir_Simulation
    properties
        ProductionSystem
        FluidModel
        DiscretizationModel
        Formulation
        TimeDriver
        Summary
    end
    methods
        function Initialize(obj)
            obj.DiscretizationModel.Initialize(obj.ProductionSystem);
            obj.ProductionSystem.Initialize(obj.DiscretizationModel, obj.FluidModel);
        end
        function RunSimulation(obj, Writer)
            obj.TimeDriver.SolveTimeDependentProblem(obj.ProductionSystem, obj.FluidModel, obj.DiscretizationModel, obj.Formulation, obj.Summary, Writer);
        end
    end
end