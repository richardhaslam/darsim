%  ADM operators handler MS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 August 2017
%Last modified: 4 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef prolongation_builder < matlab.mixin.Heterogeneous & handle
    properties
        P
        R
    end
    methods
        function obj = prolongation_builder(n)
        end
        function MsR = MsRestriction(obj, FineGrid, CoarseGrid)
            Nf = FineGrid.N;
            Nc = CoarseGrid.N;
            %% MSFV Restriction Operator
            MsR = sparse(Nc, Nf);
            for c = 1:Nc
                MsR(c, CoarseGrid.Children(c,:)) = 1;
            end
        end
    end
    methods (Abstract)
        obj = BuildStaticOperators(obj)
        obj = ADMProlongation(obj)
    end
end