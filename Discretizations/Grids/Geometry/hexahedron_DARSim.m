% Class of hexahedron for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2020-03-09
%Last modified: 2020-03-09
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
        function obj = hexahedron_DARSim(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot)
            obj.N_Vertex = 8;
            obj.N_Face = 6;
            obj.Face = tetragon_DARSim.empty;
            if nargin==8
                obj.InitializeHexahedron(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot)
            end
        end
        function InitializeHexahedron(obj,NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot)
            obj.AddCorners(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot);
            obj.InitializeFaces();
        end
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
        function InitializeFaces(obj)
            obj.Face(1) = tetragon_DARSim(obj.NW_Bot,obj.NE_Bot,obj.NE_Top,obj.NW_Top,'North');
            obj.Face(2) = tetragon_DARSim(obj.SW_Bot,obj.SE_Bot,obj.SE_Top,obj.SW_Top,'South');
            obj.Face(3) = tetragon_DARSim(obj.SE_Bot,obj.NE_Bot,obj.NE_Top,obj.SE_Top,'East');
            obj.Face(4) = tetragon_DARSim(obj.SW_Bot,obj.NW_Bot,obj.NW_Top,obj.SW_Top,'West');
            obj.Face(5) = tetragon_DARSim(obj.SW_Top,obj.SE_Top,obj.NE_Top,obj.NW_Top,'Top');
            obj.Face(6) = tetragon_DARSim(obj.SW_Bot,obj.SE_Bot,obj.NE_Bot,obj.NW_Bot,'Bottom');
        end
        function [Geostatus, IntersectPoints] = Obtain_Hexahedron_LineSegment_Intersection(obj,LineSegment)
            % The algorithm used in this section is implemented from the method found in this article:
            % http://geomalgorithms.com/a13-_intersect-4.html
            %%
            % Initializing some variables
            doNotContinue       = 0;
            Geostatus.areParallel   = NaN;
            Geostatus.areCoplanar   = NaN;
            Geostatus.haveIntersect = NaN;
            IntersectPoints         = [];
            
            dx = norm(obj.Face(1).PointM - obj.Face(2).PointM);
            dy = norm(obj.Face(3).PointM - obj.Face(4).PointM);
            dz = norm(obj.Face(5).PointM - obj.Face(6).PointM);
            Epsilon = 1e-10 * mean(dx,dy,dz);
            
            tE = 0; % for the maximum entering segment parameter
            tL = 1; % for the minimum leaving segment parameter
            for i= 1 : obj.N_Faces
                if doNotContinue
                    return;
                end
                % Checking if the normal vector is pointing outward, and if not, fixing it:
                if dot( obj.Centroid - obj.Face(i).PointM , obj.Face(i).n_vec) < 0
                    obj.Face(i).n_vec = - obj.Face(i).n_vec;
                end
                
                N = - dot( lineSegment.PointA - obj.Face(i).PointA , obj.Face(i).n_vec);
                D = dot( LineSegment.AB_vec , obj.Face(i).n_vec );
                
                if abs(D) < Epsilon
                    % Then, the lineSegment is parallel to this face
                    if N < 0
                        % Then, PointA of the lineSegment is outside of this face
                        % Leaving the entire function for this hexahedron as
                        % this lineSegment will not intersect this hexahedron.
                        Geostatus.haveIntersect = 0;
                        return;
                    else
                        
                    end
                end
                
            end
            
        end
    end
end

