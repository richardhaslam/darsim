% Class of hexahedron for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim Reservoir Simulator
% Author: Mousa HosseiniMehr
% TU Delft
% Created: 2020-03-09
% Last modified: 2020-06-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef hexahedron_DARSim < polyhedron_DARSim
    properties
        NW_Top
        SW_Top
        SE_Top
        NE_Top
        NW_Bot
        SW_Bot
        SE_Bot
        NE_Bot
        Corners
        Centroid
        Face
        Volume
    end
    methods
        %%
        function obj = hexahedron_DARSim(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot)
            obj.NumOfVertex = 8;
            obj.NumOfFace = 6;
            obj.Face = tetragon_DARSim.empty;
            if nargin >= 8
                obj.InitializeHexahedron(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot);
            end
        end
        %%
        function InitializeHexahedron(obj,NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot)
            obj.AddCorners(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot);
            obj.InitializeFaces();
        end
        %%
        function AddCorners(obj,NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot)
            obj.NW_Top = NW_Top;
            obj.SW_Top = SW_Top;
            obj.SE_Top = SE_Top;
            obj.NE_Top = NE_Top;
            obj.NW_Bot = NW_Bot;
            obj.SW_Bot = SW_Bot;
            obj.SE_Bot = SE_Bot;
            obj.NE_Bot = NE_Bot;
            obj.Corners = [NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot];
        end
        %%
        function InitializeFaces(obj)
            obj.Face(1) = tetragon_DARSim(obj.NW_Bot,obj.NE_Bot,obj.NE_Top,obj.NW_Top,'North');
            obj.Face(2) = tetragon_DARSim(obj.SW_Bot,obj.SE_Bot,obj.SE_Top,obj.SW_Top,'South');
            obj.Face(3) = tetragon_DARSim(obj.SE_Bot,obj.NE_Bot,obj.NE_Top,obj.SE_Top,'East');
            obj.Face(4) = tetragon_DARSim(obj.SW_Bot,obj.NW_Bot,obj.NW_Top,obj.SW_Top,'West');
            obj.Face(5) = tetragon_DARSim(obj.SW_Top,obj.SE_Top,obj.NE_Top,obj.NW_Top,'Top');
            obj.Face(6) = tetragon_DARSim(obj.SW_Bot,obj.SE_Bot,obj.NE_Bot,obj.NW_Bot,'Bottom');
        end
    end
end

