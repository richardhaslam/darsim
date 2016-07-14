% Output base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 8 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef output_writer < handle
    properties
        Directory
        ProblemName
        Plotter
    end
    methods
        function obj = output_writer(dir, problem)
            obj.Directory = strcat(dir, '/Output/');
            obj.ProblemName = problem;
        end
        function AddPlotter(obj, plotter)
            obj.Plotter = plotter;
        end
    end
end