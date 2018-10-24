% Run Summary for ADM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Run_Summary_ADM < Run_Summary
    properties
        ADMStats
    end
    methods
        function obj = Run_Summary_ADM(MaxNTimeSteps, couplingstats, wellsData, maxLevel)
            obj@Run_Summary(MaxNTimeSteps, couplingstats, wellsData);
            obj.ADMStats = zeros(MaxNTimeSteps, maxLevel + 2);
        end
        function SaveGridStats(obj, Ndt, DiscretizationModel)
            obj.ADMStats(Ndt,:) = [DiscretizationModel.ADMStats.N, sum(DiscretizationModel.ADMStats.N)/DiscretizationModel.ReservoirGrid.N *100];
        end
    end
end
