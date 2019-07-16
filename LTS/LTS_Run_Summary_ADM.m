%%%%
    % NOT USED ANY MORE.
%%%%%%
% Run Summary for ADM LTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% classdef LTS_Run_Summary_ADM < Run_Summary_ADM
%     properties
%         internal
%     end
%     methods
%         function obj = LTS_Run_Summary_ADM(MaxNTimeSteps, couplingstats, wellsData, maxLevel)
%             obj@Run_Summary_ADM(MaxNTimeSteps, couplingstats, wellsData, maxLevel);
%             internal =1;
%             val = 1;
%             % 2 shold be more general (time_ref problerties inside
%             % LTS_ADM_Sequatial strategy)
%             for i = 1:maxLevel
%                 val = val * 4;
%                 internal = internal + val;           
%             end
%             % Matteo: I do not fully follow what's happening here.
%             obj.internal = internal;
%             obj.CouplingStats.Complexity = zeros(MaxNTimeSteps, internal);
%         end
%     end
% end