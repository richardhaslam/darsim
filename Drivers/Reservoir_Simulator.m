% Reservoir Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
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
        end
        function BuildObjects(obj)
            obj.Builder.FindKeyWords(obj.Reader.InputMatrix, obj.Reader.SettingsMatrix);
            obj.Simulation = obj.Builder.BuildSimulation(obj.Reader.InputMatrix{1}, obj.Reader.SettingsMatrix{1});
            obj.Writer = obj.Builder.BuildWriter(obj.Simulation);
        end
        function AddSimulation(obj, simulation)
            obj.Simulation = simulation;
        end
        function AddOutput(obj, output)
            obj.Output = output;
        end
    end
end