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
        ProductionSystem
        FluidModel
        DiscretizationModel
        Formulation
        TimeDriver
        Summary
    end
    methods
        function Initialize(obj)
            obj.Formulation.SinglePhase = obj.ProductionSystem.Initialize(obj.DiscretizationModel, obj.FluidModel);
            obj.DiscretizationModel.Initialize(obj.ProductionSystem, obj.FluidModel);
            obj.ProductionSystem.InitializeWells(obj.FluidModel);
            obj.DiscretizationModel.ReservoirGrid.ComputeDepth(obj.Formulation.GravityModel.alpha);
        end
        function Run(obj, Writer)        
            disp('BEGIN TIME-DEPENDENT SIMULATION');
            disp(char(5))
            obj.TimeDriver.SolveTimeDependentProblem(obj.ProductionSystem, obj.FluidModel, obj.DiscretizationModel, obj.Formulation, obj.Summary, Writer);
            disp('end of the simulation');
        end
    end
end