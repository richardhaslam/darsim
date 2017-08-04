%  Prolongation builder MS for hyperbolic variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 August 2017
%Last modified: 4 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef prolongation_builder_MSHyperbolic < prolongation_builder
    properties

    end
    methods
        function obj = prolongation_builder_MSHyperbolic(n, cf)
            obj@prolongation_builder(n)
            obj.P = cell(1, n);
        end
        function ADMProl = ADMProlongation(obj, ADMGrid, FineGrid, CoarseGrid, ADMRest)
            
        end
    end
end