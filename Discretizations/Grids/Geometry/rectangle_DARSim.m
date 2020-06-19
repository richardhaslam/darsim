% class rectangle for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim 2 Reservoir Simulator
% Author: Mousa HosseiniMehr
% TU Delft
% Created: 2020-06-17
% Last modified: 2020-06-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef rectangle_DARSim < tetragon_DARSim
    properties
    end
    methods
        %%
        function obj = rectangle_DARSim(pointA,pointB,pointC,pointD,position)
            obj.NumOfVertex = 4;
            if nargin >= 4
                obj.InitializeTetragon(pointA,pointB,pointC,pointD);
            end
            if nargin >= 5
                obj.Position = position;
            end
            % We need to create a function that checks if this is a rectangle or not.
        end
    end
end