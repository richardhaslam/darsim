%  FS discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FS_Discretization_model < Discretization_model
    properties
        
    end
    methods
        function InitializeMapping(obj, ProductionSystem, FluidModel)
            % virtual call
        end
        function AverageSolOnCoarseBlocks(obj, Status, FluidModel, Formulation)
            % virtual call
        end
        function SelectADMGrid(obj, ProductionSystem)
            % virtual call
        end
        function AddWellsToInitialPressure(obj, ProductionSystem, FluidModel)
            % virtual call
        end
    end
end