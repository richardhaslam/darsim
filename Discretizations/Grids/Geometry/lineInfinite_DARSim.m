% Class of infinite line for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim Reservoir Simulator
% Author: Mousa HosseiniMehr
% TU Delft
% Created: 2020-03-09
% Last modified: 2020-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef lineInfinite_DARSim < handle
    properties
        PointA
        unitVec
        Eq       % The Eq of line is applied as follows x = a*t+x0; y = b*t+y0; z = c*t+z0;
        isParallelToXYZAxis
        isParallelToXYZSurface
    end
    methods
        %%
        function obj = lineInfinite_DARSim(unitVec,pointA)
            if nargin >= 1
                obj.InitializeLineInfinite(unitVec);
            end
            if nargin >= 2
                obj.AddPoint(pointA);
            end
        end
        %%
        function InitializeLineInfinite(obj,unitVec)
            obj.unitVec = unitVec / norm(unitVec);
        end
        %%
        function AddPoint(obj,pointA)
            obj.PointA = pointA;
        end
        %%
        function AddEquation(obj,pointA)
            if nargin==2 && isempty(obj.PointA)
                obj.AddPoint(pointA);
            end
            obj.Eq.a  = obj.unitVec(1);
            obj.Eq.b  = obj.unitVec(2);
            obj.Eq.c  = obj.unitVec(3);
            obj.Eq.x0 = obj.PointA(1);
            obj.Eq.y0 = obj.PointA(2);
            obj.Eq.z0 = obj.PointA(3);
        end
        %%
        function [Geostatus, intersectPoint] = Obtain_LineInfinite_LineInfinite_Intersection(obj, Line, Epsilon)
            Geostatus.areParallel = NaN;
            Geostatus.areCollinear = NaN;
            Geostatus.haveIntersect = NaN;
            intersectPoint = [];
            
            % Obtaining the geostatus between the lines
            if ( norm( cross(obj.unitVec, Line.unitVec) ) < Epsilon )
                % The lines are either parallel or even collinear.
                Geostatus.areParallel = 1;
                
                % checking if there is any distance between the prarallel line segment
                parallelDistance = norm( cross( (Line.PointA - obj.PointA) , obj.unitVec ) );
                if parallelDistance > Epsilon
                    % The lines are parallel, no intersection occurs.
                    Geostatus.areCollinear = 0;
                    Geostatus.haveIntersect = 0;
                else
                    % The lines are collinear.
                    Geostatus.areCollinear = 1;
                end
                
            else
                % The lines are either intersected or skew.
                Geostatus.areParallel  = 0;
                Geostatus.areCollinear = 0;
                L1A = obj.PointA  - obj.unitVec  * 1e3;
                L1B = obj.PointA  + obj.unitVec  * 1e3;
                L2A = Line.PointA - Line.unitVec * 1e3;
                L2B = Line.PointA + Line.unitVec * 1e3;
                V1 = L1B - L1A;
                V2 = L2B - L2A;
                t = det( [ (L2A-L1A)', V2' , cross(V1,V2)' ] ) / norm(cross(V1,V2))^2;
                s = det( [ (L2A-L1A)', V1' , cross(V1,V2)' ] ) / norm(cross(V1,V2))^2;
                
                if norm( (L1A + t*V1) - (L2A + s*V2) ) < Epsilon
                    % Intersection occurs
                    Geostatus.haveIntersect  = 1;
                    intersectPoint = L1A + V1.*t;
                else
                    % The lines are skew
                    Geostatus.haveIntersect  = 0;
                end
            end
        end
        %%
        function IsPointOnTheInfiniteLine = Is_Point_On_The_InfiniteLine(obj, Point, Epsilon)
            % The point P lies on an infinite line (with P0 being a point defined on that line)
            % only if the angle between the unit vector of the line and PP0 is zero or 180 degree.
            % This means that if tge cross product of these two vectors is zero, the point is on the line.
            if cross( Point - obj.PointA , obj.unitVec ) < Epsilon
                IsPointOnTheInfiniteLine = 1;
            else
                IsPointOnTheInfiniteLine = 0;
            end
        end
        %%
        function Distance = Obtain_Distance_Point_from_InfiniteLine(obj, Points)
            Distance = vecnorm( cross(Points-obj.PointA , obj.unitVec.*ones(size(Points))) , 2 , 2 ) ./ vecnorm( obj.unitVec , 2 , 2 );
        end
    end
end