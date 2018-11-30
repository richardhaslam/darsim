% Run Summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Run_Summary < handle
    properties
        WellsData
        CouplingStats
        NumberTimeSteps
        Time
        NofPhases
        NofComponents
    end
    methods
        function obj = Run_Summary(MaxNTimeSteps, couplingstats, wellsData)
            obj.WellsData = wellsData;
            obj.CouplingStats = couplingstats;
            obj.Time = zeros(MaxNTimeSteps, 1);
        end
        function SaveWellsData (obj, Ndt, Inj, Prod, dt)
            obj.WellsData.UpdateInjectionCurve(Ndt, Inj, dt);
            obj.WellsData.UpdateProductionCurve(Ndt, Prod, dt);
        end 
        function SaveGridStats(obj, Ndt, DiscretizationModel)
            % virtual call
            
        end
    end
end
