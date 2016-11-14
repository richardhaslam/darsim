%  FS discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FS_Discretization_model < Discretization_model
    properties
        
    end
    methods
        function obj = FS_Discretization_model(nx, ny, nz)
            obj@Discretization_model(nx, ny, nz);
        end
        function InitializeMapping(obj, ProductionSystem, FluidModel)
            % virtual call
        end
        function AverageSolOnCoarseBlocks(obj, Status, FluidModel, Formulation)
            % virtual call
        end
    end
end