%  FS discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 21 September 2017
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
    end
end