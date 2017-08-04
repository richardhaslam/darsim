%  Prolongation builder constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 August 2017
%Last modified: 4 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef prolongation_builder_constant < prolongation_builder
    properties
    end
    methods
        function obj = prolongation_builder_constant(n)
            obj@prolongation_builder(n)
        end
        function BuildStaticOperators(obj, CoarseGrid, FineGrid, maxLevel, K, S, FluidModel)
            % Build Restriction and Prolongation operators for static grids
        end
        function ADMProl = ADMProlongation(obj, ADMGrid, FineGrid, CoarseGridid, ADMRest)
            % Since it s constant interpolation it is just transpose(R)
            ADMProl = ADMRest';
        end
    end
end