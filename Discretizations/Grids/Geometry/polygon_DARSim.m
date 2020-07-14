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
        Vertex
        Centroid
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
        function AddVertices(obj,points)
            if size(points,1) ~= obj.NumOfVertex
                error('ADMIRE Error: the number of points given does not match the number of vertices');
            end
            obj.Vertex = points;
        end
        %%
        function OrderCornersClockwise(obj)
            % ordering the corners clockwise
            % (https://nl.mathworks.com/matlabcentral/answers/429265-how-to-order-vertices-of-a-flat-convex-polygon-in-3d-space-along-the-edge)
            xyz = obj.Vertex;
            xyzc = mean(xyz,1);
            P = xyz - xyzc;
            [~,~,V] = svd(P,0);
            [~,is] = sort(atan2(P*V(:,1),P*V(:,2)));
            obj.Vertex = obj.Vertex(is,:);
        end
        %%
        function ObtainCentroid(obj)
            obj.Centroid = mean(obj.Vertex,1);
        end
        %%
        function isInside = Is_Point_Inside_Polygin(obj, point , Epsilon)
            Line1 = lineSegment_DARSim(obj.Centroid , point);
            for n = 1 : obj.NumOfVertex
                Line2 = lineSegment_DARSim( obj.Vertex(n,:) , obj.Vertex(mod((end+n)-1,end)+1,:) );
                [Geostatus, IntersectPoint] = Line1.Obtain_LineSegment_LineSegment_Intersection( Line2, Epsilon );
                if     ( Geostatus.haveIntersect == 0 )             ,  isInside = 1;
                elseif ( norm (IntersectPoint - point ) < Epsilon ) ,  isInside = 1;
                else                                                ,  isInside = 0;  return;  end    
            end
        end
        %%
        function [Geostatus, IntersectPoint] = Obtain_Polygon_LineSegment_Intersection(obj, LineSegment, Epsilon)
            % 1. First, we check if the plane of this polygon intersect with the line segment
            [Geostatus, IntersectPoint] = obj.Obtain_PlaneInfinite_LineSegment_Intersection(LineSegment, Epsilon);
            
            if Geostatus.haveIntersect == 1 && ~isempty(IntersectPoint)
                % 2. Now, we check if the intersection point lies inside the polygon
                PointIsInsidePolygon = obj.Is_Point_Inside_Polygin(IntersectPoint , Epsilon);
                
                if PointIsInsideLineSegment && PointIsInsidePolygon
                    % The information is already correctly reported with [Geostatus, intersectPoint].
                    % No further action is needed.
                else
                    Geostatus.haveIntersect = NaN;
                    IntersectPoint = [];
                end
            end
        end
    end
end