% Class of polygon for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim Reservoir Simulator
% Author: Mousa HosseiniMehr
% TU Delft
% Created: 2020-03-09
% Last modified: 2020-06-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef polygon_DARSim < planeInfinite_DARSim
    properties
        NumOfVertex
        Position
    end
    methods
        %%
        function obj = polygon_DARSim(numOfVertex,n_vec,point)
            if nargin >= 1
                obj.NumOfVertex = numOfVertex;
            end
            if nargin >= 2
                obj.InitializePolygon(n_vec,point);
            end
        end
        %%
        function InitializePolygon(obj,n_vec,point)
            obj.InitializePlaneInfinite(n_vec,point);
        end
        %%
    end
end