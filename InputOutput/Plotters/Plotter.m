% Plotter base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Plotter < handle
    properties
        VTKindex
    end
    methods (Abstract)
       obj = PlotSolution(obj) 
       obj = PlotPermeability(obj)
       obj = PlotWells(obj)
       obj = PlotADMGrid(obj)
    end
    methods
%         function PlotInitialStatus(obj, ProductionSystem, DiscretizationModel)
%             % Plot reservoir Solution
%             obj.PlotSolution(ProductionSystem.Reservoir.State, DiscretizationModel.ReservoirGrid);
%             obj.PlotPermeability(DiscretizationModel.ReservoirGrid, ProductionSystem.Reservoir.K);
%             % Plot Fractures
%             
%             
%         end
    end
end