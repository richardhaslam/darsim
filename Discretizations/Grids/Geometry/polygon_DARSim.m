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
        function isInside = Is_Point_Inside_Polygon(obj, point , Epsilon)
            Line1 = lineSegment_DARSim(obj.Centroid , point);
            for n = 1 : obj.NumOfVertex
                Line2 = lineSegment_DARSim( obj.Vertex(n,:) , obj.Vertex(mod((end+n)-1,end)+1,:) );
                [Geostatus, IntersectPoint] = Line1.Obtain_LineSegment_LineSegment_Intersection( Line2, Epsilon );
                if     ( Geostatus.haveIntersect == 0 )                                         ,  isInside = 1;
                elseif ~isempty(IntersectPoint) && ( norm (IntersectPoint - point ) < Epsilon ) ,  isInside = 1;
                else                                                                            ,  isInside = 0;  return;  end
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
        %%
        function [Geostatus, IntersectPoints] = Obtain_Polygon_Polygon_Intersection(obj, Polygon, Epsilon )
            Geostatus.areParallel   = NaN;
            Geostatus.areCoplanar   = NaN;
            Geostatus.haveIntersect = NaN;
            IntersectPoints         = [];
            
            % Obtaining the geostatus between the polygons
            if norm( cross( obj.nVec , Polygon.nVec ) ) < Epsilon
                % The polygons are parallel
                Geostatus.areParallel = 1;
                % If the equations of both polygons are multiple of each
                % other, then they are coplanar:
                % To check this, we put a point from the 1st polygon into
                % 2nd polygon's equation. If the equation holds, it
                % means that the point lies within both polygons resulting in
                % coplanarity of these two polygons.
                % if a2*x1 + b2*y1 + c2*z1 = d2, then the polygons are coplanar.
                if abs( Polygon.Eq.a * obj.Centroid(1) + ...
                        Polygon.Eq.b * obj.Centroid(2) + ...
                        Polygon.Eq.c * obj.Centroid(3) - ...
                        Polygon.Eq.d                     ) < Epsilon
                    % The polygons are coplanar
                    Geostatus.areCoplanar = 1;
                    
                    % Checking the intersections between each two edges of polygons
                    for i = [1,2,3,4]
                        intersectNr = 0;
                        for j = [2,3,4,1]
                            Line1 = lineSegment_DARSim(obj.Vertex(i,:)     ,obj.Vertex(j,:)     );
                            Line2 = lineSegment_DARSim(Polygon.Vertex(i,:),Polygon.Vertex(j,:));
                            [lineGeostatus, lineIntersectPoint] = Line1.Obtain_LineSegment_LineSegment_Intersection(Line2, Epsilon);
                            if lineGeostatus.haveIntersect == 1
                                Geostatus.haveIntersect = 1;
                                intersectNr = intersectNr + 1;
                                IntersectPoints = [IntersectPoints ; lineIntersectPoint];
                            end
                            if intersectNr == 2,  continue;  end  % each edge of 1st plane segment can have intersection with maximum two edges of the 2nd plane segment
                        end
                    end
                    
                else
                    % The polygons are not coplanar
                    Geostatus.areCoplanar = 0;
                    Geostatus.haveIntersect = 0;
                end
                
            else
                % The polygons are not parallel and have intersection line
                Geostatus.areParallel = 0;
                Geostatus.areCoplanar = 0;
                Geostatus.haveIntersect = 1;
                
                % Obtaining the unit vector of the intersection line
                intLine_unitVec = cross( obj.nVec , Polygon.nVec );
                intLine_unitVec = intLine_unitVec / norm(intLine_unitVec);
                
                % Obtaining a Point on the intersection line with Z=0 if possible
                P2M = mean(Polygon.Vertex);
                Axy = [ obj.nVec(1) , obj.nVec(2) ; Polygon.nVec(1) , Polygon.nVec(2) ];
                RHS = [ obj.Eq.d      - obj.nVec(3)*P2M(3)
                        Polygon.Eq.d - Polygon.nVec(3)*P2M(3) ];
                if abs(det(Axy)) > Epsilon
                    intLine_Point0    = zeros(1,3);
                    Unknowns          = Axy \ RHS;
                    intLine_Point0(1) = Unknowns(1);
                    intLine_Point0(2) = Unknowns(2);
                    intLine_Point0(3) = P2M(3);
                else
                    % Or, obtaining a Point on the intersection line with Y=0 if possible
                    Axz = [ obj.nVec(1) , obj.nVec(3) ; Polygon.nVec(1) , Polygon.nVec(3) ];
                    RHS = [ obj.Eq.d      - obj.nVec(2)*P2M(2)
                            Polygon.Eq.d - Polygon.nVec(2)*P2M(2) ];
                    if abs(det(Axz)) > Epsilon
                        intLine_Point0    = zeros(1,3);
                        Unknowns          = Axz \ RHS;
                        intLine_Point0(1) = Unknowns(1);
                        intLine_Point0(3) = Unknowns(2);
                        intLine_Point0(2) = P2M(2);
                    else
                        % Or, obtaining a Point on the intersection line with X=0 if possible
                        Ayz = [ obj.nVec(2) , obj.nVec(3) ; Polygon.nVec(2) , Polygon.nVec(3) ];
                        RHS = [ obj.Eq.d      - obj.nVec(1)*P2M(1)
                                Polygon.Eq.d - Polygon.nVec(1)*P2M(1) ];
                        intLine_Point0    = zeros(1,3);
                        Unknowns       = Ayz \ RHS;
                        intLine_Point0(2) = Unknowns(1);
                        intLine_Point0(3) = Unknowns(2);
                        intLine_Point0(1) = P2M(1);
                    end
                end
                
                % Assuming two end-Points for the intersection line
                Radius1 = norm( mean(    obj.Vertex) -     obj.Vertex(1,:) );
                Radius2 = norm( mean(Polygon.Vertex) - Polygon.Vertex(1,:) );
                RelativeLength = max( Radius1, Radius2 );
                intLine_A = intLine_Point0 - intLine_unitVec * RelativeLength * 1e5;
                intLine_B = intLine_Point0 + intLine_unitVec * RelativeLength * 1e5;
                intersectionLine = lineSegment_DARSim(intLine_A ,intLine_B);
                
                % Checking the intersection between the obtained line segment and each side of each polygon
                AllPoints = [obj.Vertex ; Polygon.Vertex];
                ind1 = [  1 : obj.NumOfVertex+Polygon.NumOfVertex  ];
                ind2 = [  2 : obj.NumOfVertex , 1  ,  obj.NumOfVertex+2:obj.NumOfVertex+Polygon.NumOfVertex , obj.NumOfVertex+1  ];
                intersectNr = 0;
                for i = 1 : obj.NumOfVertex+Polygon.NumOfVertex
                    Line_temp = lineSegment_DARSim( AllPoints(ind1(i),:) , AllPoints(ind2(i),:) );
                    [lineGeostatus, lineIntersectPoint] = Line_temp.Obtain_LineSegment_LineSegment_Intersection( intersectionLine, Epsilon );
                    if lineGeostatus.haveIntersect == 1
                        Geostatus.haveIntersect = 1;
                        intersectNr = intersectNr + 1;
                        IntersectPoints = [IntersectPoints ; lineIntersectPoint];
                    end
                    if intersectNr == 4,  continue;  end % intersection line can have intersection with maximum two edges of each plane segment (four intersections max)
                end
                
            end
            
        end
    end
end