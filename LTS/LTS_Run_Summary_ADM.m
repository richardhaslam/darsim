% Run Summary for ADM LTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef LTS_Run_Summary_ADM < Run_Summary_ADM
    properties
        LTSIter
        internal
    end
    methods
        function obj = LTS_Run_Summary_ADM(MaxNTimeSteps, couplingstats, wellsData, maxLevel)
            obj@Run_Summary_ADM(MaxNTimeSteps, couplingstats, wellsData, maxLevel);
            internal = 0;
            val = 1;
            % 2 shold be more general (time_ref problerties inside
            % LTS_ADM_Sequatial strategy)
            for i = 1:maxLevel
                val = val * 2;
                internal = internal + val;           
            end
            obj.LTSIter = zeros(MaxNTimeSteps, internal);
            obj.internal = internal;
        end
        function SaveLTSiter(obj, Ndt, LTS_iters)
            LTS_iters2 = zeros(1,obj.internal);
            LTS_iters2(1,1:size(LTS_iters,2)) = LTS_iters;
            obj.LTSIter(Ndt-1,:) = LTS_iters2;
        end
    end
end