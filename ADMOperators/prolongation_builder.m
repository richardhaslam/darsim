%  ADM operators handler MS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
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
            %% MSFV Restriction Operator for F-AMS
            MsR = [];
            for i=1:length(FineGrid)
                Nf = FineGrid(i).N;
                Nc = CoarseGrid(i).N;
                if FineGrid(i).CoarseLevel > 0 && FineGrid(i).hasCoarseNodes == 0
                    Nf = 0;
                end
                if CoarseGrid(i).hasCoarseNodes == 0
                    Nc = 0;
                end
                %% MSFV Restriction Operator
                R = sparse([],[],[],Nc, Nf, Nf);
                for c = 1:Nc
                    R(c, CoarseGrid(i).Children{c,1}) = 1;
                end
                MsR = blkdiag(MsR, R);
                
            end
        end
    end
    methods (Abstract)
        obj = BuildStaticOperators(obj)
        obj = ADMProlongation(obj)
    end
end