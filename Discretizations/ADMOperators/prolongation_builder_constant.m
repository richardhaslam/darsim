%  Prolongation builder constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 August 2017
%Last modified: 7 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef prolongation_builder_constant < prolongation_builder
    properties
    end
    methods
        function obj = prolongation_builder_constant(n)
            obj@prolongation_builder(n)
        end
        function BuildStaticOperators(ProductionSystem, FluidModel, ReservoirGrid, FracturesGrid, CrossConnections, maxLevel, CoarseGrid)
            % Build Restriction and Prolongation operators for static grids
        end
        function UpdateProlongationOperator(obj, FineGrid, CoarseGrid, ProductionSystem)
            % For constant bf no update is necessary
        end
        function ADMProl = ADMProlongation(obj, ADMGrid, GlobalGrids, ADMRest)
            % Since it s constant interpolation it is just transpose(R)
            ADMProl = ADMRest';
        end
        function AverageMassOnCoarseBlocks(obj, Formulation, ProductionSystem, FineGrid, FluidModel, ADMRest)
            % Average mass on a coarse block
            Formulation.AverageMassOnCoarseBlocks(ProductionSystem, FineGrid, FluidModel, ADMRest);
        end
    end
end