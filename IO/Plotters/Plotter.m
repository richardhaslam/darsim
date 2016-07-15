% Plotter base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 11 July 2016
%Last modified: 11 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Plotter < handle
    properties
    end
    methods (Abstract)
       obj = PlotSolution(obj) 
       obj = PlotPermeability(obj)
       obj = PlotWells(obj)
       obj = PlotADMGrid(obj)
    end
    methods
        function PlotInitialStatus(obj, Status, Grid, K, Inj, Prod)
            % Plot reservoir Solution
            obj.PlotSolution(Status, Grid);
            obj.PlotPermeability(Grid, K);
            % Plot Fractures
            
            % Plot Wells
            obj.PlotWells(Inj, Prod);
        end
    end
end