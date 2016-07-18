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
        function PlotInitialStatus(obj, ProductionSystem, DiscretizationModel)
            % Plot reservoir Solution
            obj.PlotSolution(ProductionSystem.Reservoir.State, DiscretizationModel.ReservoirGrid);
            obj.PlotPermeability(DiscretizationModel.ReservoirGrid, ProductionSystem.Reservoir.K);
            % Plot Fractures
            
            % Plot Wells
            obj.PlotWells(ProductionSystem.Wells.Inj, ProductionSystem.Wells.Prod, DiscretizationModel.ReservoirGrid);
        end
    end
end