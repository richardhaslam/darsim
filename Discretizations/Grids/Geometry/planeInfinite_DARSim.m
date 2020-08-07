% Class of infinite plane for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim Reservoir Simulator
% Author: Mousa HosseiniMehr
% TU Delft
% Created: 2020-03-09
% Last modified: 2020-06-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef planeInfinite_DARSim < handle
    properties
        nVec
        Eq
        isParallelToReferenceSurface
    end
    methods
        %%
        function obj = planeInfinite_DARSim(n_vec,point)
            if nargin >= 1
                obj.nVec = n_vec / norm(n_vec);
                obj.nVec( isnan(obj.nVec) ) = 0;
            end
            if nargin >= 2
                obj.InitializePlaneInfinite(n_vec,point);
            end
        end
        %%
        function InitializePlaneInfinite(obj,n_vec,point)
            % Eq of the plane follows ax+by+cz=d.
            obj.nVec = n_vec / norm(n_vec);
            obj.nVec( isnan(obj.nVec) ) = 0;
            obj.Eq.a = obj.nVec(1);
            obj.Eq.b = obj.nVec(2);
            obj.Eq.c = obj.nVec(3);
            obj.Eq.d = obj.Eq.a * point(1) + obj.Eq.b * point(2) + obj.Eq.c * point(3);
        end
        %%
        function CheckIfParallelToReferenceSurface(obj,Epsilon)
            if     ( abs(obj.Eq.b) < Epsilon ) && ( abs(obj.Eq.c) < Epsilon ), obj.isParallelToReferenceSurface = "AlongYZ";
            elseif ( abs(obj.Eq.a) < Epsilon ) && ( abs(obj.Eq.c) < Epsilon ), obj.isParallelToReferenceSurface = "AlongXZ";
            elseif ( abs(obj.Eq.a) < Epsilon ) && ( abs(obj.Eq.b) < Epsilon ), obj.isParallelToReferenceSurface = "AlongXY";
            else                                                             , obj.isParallelToReferenceSurface = "None";
            end
        end
        %%
        function [Geostatus, IntersectPoint] = Obtain_PlaneInfinite_LineInfinite_Intersection(obj, Line, Epsilon)
            Geostatus.areParallel = NaN;
            Geostatus.areCoplanar = NaN;
            Geostatus.haveIntersect = NaN;
            IntersectPoint = [];
            
            % Obtaining the geostatus between the line and the plane
            if abs( dot(obj.nVec,Line.unitVec) ) < Epsilon
                % The line and the plane are parallel
                Geostatus.areParallel = 1;
                if abs( obj.Eq.a*Line.Eq.x0 + obj.Eq.b*Line.Eq.y0 + obj.Eq.c*Line.Eq.z0 - obj.Eq.d ) < Epsilon
                    % The line and the plane are coplanar
                    Geostatus.areCoplanar = 1;
                else
                    % The line and the plane are not coplanar
                    Geostatus.areCoplanar = 0;
                    Geostatus.haveIntersect = 0;
                end
                
            else
                % The line and the plane are not parallel and have
                % intersection
                Geostatus.areParallel = 0;
                Geostatus.areCoplanar = 0;
                Geostatus.haveIntersect = 1;
                
                % Line Eq is: x=a't+x0, y=b't+y0, z=c't+z0.
                % Plane Eq is: ax+by+cz=d.
                % By substituting line Eq into plane Eq, the
                % parameter "t" is obtained:
                % a(a't+x0)+b(b't+y0)+c(c't+z0)=d
                % ===> t = ( d-ax0-by0-cz0 ) / ( a.a' + b.b' + c.c' );
                % Afterwards, by putting the value of "t" parameter in the
                % line Eq, the [x,y,z] of the intersection point is
                % obtained.
                
                t = ( obj.Eq.d - obj.Eq.a*Line.Eq.x0 - obj.Eq.b*Line.Eq.y0 - obj.Eq.c*Line.Eq.z0 ) ...
                             / ( obj.Eq.a*Line.Eq.a' + obj.Eq.b*Line.Eq.b' + obj.Eq.c*Line.Eq.c' );
                
                x = Line.Eq.a * t + Line.Eq.x0;
                y = Line.Eq.b * t + Line.Eq.y0;
                z = Line.Eq.c * t + Line.Eq.z0;
                IntersectPoint = [x,y,z];
            end
        end
        %%
        function [Geostatus, IntersectPoint] = Obtain_PlaneInfinite_LineSegment_Intersection(obj, LineSegment, Epsilon)
            % 1. First, we check if the plane of this polygon intersect with the infinite extenstion of the line segment
            [Geostatus, IntersectPoint] = obj.Obtain_PlaneInfinite_LineInfinite_Intersection(LineSegment, Epsilon);
            if size(IntersectPoint,1) > 1
                error('ADMIRE Error: The intersection between a plane and a line cannot be more than one point.');
            end
            
            if Geostatus.haveIntersect == 1 && ~isempty(IntersectPoint)
                % 2. Now, we check if the intersection point lies inside the line segment
                % Point C is inside the line segment AB only if dot product of AC and BC is negative.
                if dot( IntersectPoint-LineSegment.PointA , IntersectPoint-LineSegment.PointB ) <= 0
                    PointIsInsideLineSegment = 1;
                else
                    PointIsInsideLineSegment = 0;
                end
                
                if PointIsInsideLineSegment
                    % The information is already correctly reported with [Geostatus, intersectPoint].
                    % No further action is needed.
                else
                    Geostatus.haveIntersect = NaN;
                    IntersectPoint = [];
                end
            end
        end
        %%
        function [Geostatus, intersectLine] = Obtian_PlaneInfinite_PlaneInfinite_Intersection(obj, Plane, Epsilon)
            Geostatus.areParallel = NaN;
            Geostatus.areCoplanar = NaN;
            Geostatus.haveIntersect = NaN;
            intersectLine = [];
            % Obtaining the geostatus between the planes
            if norm( cross( obj.nVec , Plane.nVec ) ) < Epsilon
                % The planes are parallel
                Geostatus.areParallel = 1;
                % If the equations of both planes are multiple of each
                % other, then they are coplanar:
                % d2 * (a1/a2) - d1 = 0 or d2 * (b1/b2) - d1 = 0 or d2 * (c1/c2) - d1 = 0
                if abs( Plane.Eq.d * (obj.Eq.a/Plane.Eq.a) - obj.Eq.d ) < Epsilon || ...
                   abs( Plane.Eq.d * (obj.Eq.b/Plane.Eq.b) - obj.Eq.d ) < Epsilon || ...
                   abs( Plane.Eq.d * (obj.Eq.c/Plane.Eq.c) - obj.Eq.d ) < Epsilon
                    % The planes are coplanar
                    Geostatus.areCoplanar = 1;
                else
                    % The planes are not coplanar
                    Geostatus.areCoplanar = 0;
                    Geostatus.haveIntersect = 0;
                end
            else
                % The planes are not parallel and have intersection line
                Geostatus.areParallel = 0;
                Geostatus.areCoplanar = 0;
                Geostatus.haveIntersect = 1;
                unit_vec = cross(obj.nVec,Plane.nVec);
                unit_vec = unit_vec / norm(unit_vec);
                
                % Obtaining a Point on the intersection line:
                Axy = [ obj.Eq.a , obj.Eq.b ; Plane.Eq.a , Plane.Eq.b ];
                Axz = [ obj.Eq.a , obj.Eq.c ; Plane.Eq.a , Plane.Eq.c ];
                Ayz = [ obj.Eq.b , obj.Eq.c ; Plane.Eq.b , Plane.Eq.c ];
                RHS = [ obj.Eq.d ; Plane.Eq.d ];
                % with Z=0 if possible
                if abs(det(Axy)) > Epsilon
                    xy = Axz \ RHS;
                    intersectPoint = [xy(1);xy(2);0];
                % Obtaining a Point on the intersection line with Y=0 if possible    
                elseif abs(det(Axz)) > Epsilon
                    xz = Axz \ RHS;
                    intersectPoint = [xz(1);0;xz(2)];
                else % abs(det(Ayz)) > Epsilon
                    yz = Ayz \ RHS;
                    intersectPoint = [0;yz(1);yz(2)];
                end
                intersectLine = lineInfinite_DARSim(unit_vec,intersectPoint);   
            end
        end
    end
end