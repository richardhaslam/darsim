% Class of polygon for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2020-03-09
%Last modified: 2020-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef polygon_DARSim < plane_DARSim
    properties
        N_Vertex
        Pos
    end
    methods
        function obj = polygon_DARSim(N_vertex,n_vec,point,pos)
            %%
            if nargin>=1
                obj.N_Vertex = N_vertex;
            end
            if nargin>=3
                obj.InitializePolygon(n_vec,point);
            end
            if nargin==4
                obj.AddPosition(pos);
            end
        end
        function InitializePolygon(obj,n_vec,point)
            obj.InitializePlane(n_vec,point);
        end
        function AddPosition(obj,pos)
            obj.Pos = pos;
        end
    end
end