% Simulation class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 17 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Reservoir_Simulation
    properties
        Initializer
        ProductionSystem
        FluidModel
        DiscretizationModel
        Formulation
        TimeDriver
        Summary
    end
    methods
        function Initialize(obj)
            obj.DiscretizationModel.Initialize(obj.ProductionSystem, obj.Formulation);
            obj.Initializer.InitializeProductionSystem(obj.ProductionSystem, obj.FluidModel, obj.Formulation, obj.DiscretizationModel);
            obj.ProductionSystem.InitializeWells(obj.FluidModel, obj.Formulation.GravityModel, obj.DiscretizationModel);
            obj.DiscretizationModel.InitializeMapping(obj.ProductionSystem, obj.FluidModel);
            % For multiscale to avoid having well functions
            % obj.DiscretizationModel.AddWellsToInitialPressure(obj.ProductionSystem, obj.FluidModel);
        end
        function Run(obj, Writer)
            disp('BEGIN TIME-DEPENDENT SIMULATION');
            disp(char(5))
            obj.TimeDriver.SolveTimeDependentProblem(obj.ProductionSystem, obj.FluidModel, obj.DiscretizationModel, obj.Formulation, obj.Summary, Writer);
            disp('end of the simulation');
        end
    end
end