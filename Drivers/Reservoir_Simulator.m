% Reservoir Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 18 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Reservoir_Simulator < handle
    properties
        Reader
        Builder
        Simulation
        Writer
    end
    methods
        function obj = Reservoir_Simulator(Directory, File)
            obj.Reader = reader(Directory, File);
            obj.Builder = builder();
            obj.Simulation = Reservoir_Simulation();
        end
        function BuildObjects(obj)
            obj.Builder.FindKeyWords(obj.Reader.InputMatrix, obj.Reader.SettingsMatrix);
            obj.Simulation = obj.Builder.BuildSimulation(obj.Reader.InputMatrix{1}, obj.Reader.SettingsMatrix{1});
            obj.Writer = obj.Builder.BuildWriter(obj.Reader.Directory, obj.Simulation); 
        end
        function Run(obj)
            % Initialize Simulation
            obj.Simulation.Initialize();
            % Plot initial state of the reservoir
            obj.Writer.Plotter.PlotInitialStatus(obj.Simulation.ProductionSystem, obj.Simulation.DiscretizationModel);
            % Run simulation
            obj.Simulation.Run(obj.Writer);
        end
        function OutputResults(obj)
            obj.Writer.WriteSummary(obj.Simulation.Summary);
        end
    end
end