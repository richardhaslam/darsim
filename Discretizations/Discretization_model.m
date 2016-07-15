%  Discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Discretization_model < handle
    properties
        ReservoirGrid
        FracturesGrid
    end
    methods
        function obj = Discretization_model(nx, ny, nz)
            obj.ReservoirGrid = grid_cartesian(nx, ny, nz);
        end
        function Initialize(obj, ProductionSystem)
            obj.ReservoirGrid.Initialize(ProductionSystem.Reservoir);
        end
    end
end