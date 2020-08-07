% Class of line segment for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim Reservoir Simulator
% Author: Mousa HosseiniMehr
% TU Delft
% Created: 2020-03-09
% Last modified: 2020-06-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef lineSegment_DARSim < lineInfinite_DARSim
    properties
        PointB
        PointM
        Points
        AB_vec
    end
    methods
        %%
        function obj = lineSegment_DARSim(pointA,pointB)
            obj.InitializeLineInfinite(pointB-pointA);
            obj.PointA = pointA;
            obj.PointB = pointB;
            obj.PointM = (pointA+pointB)/2;
            obj.AB_vec = pointB - pointA;
        end
        %%
        function [Geostatus, intersectPoint] = Obtain_LineSegment_LineSegment_Intersection(obj, Line, Epsilon)
            Geostatus.areParallel = NaN;
            Geostatus.areCollinear = NaN;
            Geostatus.areSkew = NaN;
            Geostatus.haveIntersect = NaN;
            intersectPoint = [];
            
            L1A = obj.PointA; L2A = Line.PointA;
            L1B = obj.PointB; L2B = Line.PointB;
            V1 = L1B - L1A;
            V2 = L2B - L2A;
            
            % Obtaining the geostatus between the lines
            if ( norm( cross(obj.unitVec, Line.unitVec) ) < Epsilon )
                % The lines are either parallel or even collinear.
                Geostatus.areParallel = 1;
                Geostatus.areSkew = 0;
                
                % checking if there is any distance between the prarallel line segment
                parallelDistance = norm( cross( (Line.PointA - obj.PointA) , obj.unitVec ) );
                if parallelDistance > Epsilon
                    % The lines are parallel, no intersection occurs.
                    Geostatus.areCollinear = 0;
                    Geostatus.haveIntersect = 0;
                else
                    % The lines are collinear.
                    Geostatus.areCollinear = 1;
                    intersectNr = 0;

                    % Checking if point L2A is an intersection or not
                    arcCos_Theta_L2A_1 = dot( L1B-L1A , L2A-L1A ) / ( norm(L1B-L1A) * norm(L2A-L1A) );
                    arcCos_Theta_L2A_2 = dot( L1A-L1B , L2A-L1B ) / ( norm(L1A-L1B) * norm(L2A-L1B) );
                    if ( arcCos_Theta_L2A_1 * arcCos_Theta_L2A_2 > Epsilon ) && ( intersectNr <= 2 )
                        intersectNr = intersectNr + 1;
                        Geostatus.haveIntersect = 1;
                        intersectPoint = [ intersectPoint ; L2A ];
                    end
                    
                    % Checking if point L2B is an intersection or not
                    arcCos_Theta_L2B_1 = dot( L1B-L1A , L2B-L1A ) / ( norm(L1B-L1A) * norm(L2B-L1A) );
                    arcCos_Theta_L2B_2 = dot( L1A-L1B , L2B-L1B ) / ( norm(L1A-L1B) * norm(L2B-L1B) );
                    if ( arcCos_Theta_L2B_1 * arcCos_Theta_L2B_2 > Epsilon ) && ( intersectNr <= 2 )
                        intersectNr = intersectNr + 1;
                        Geostatus.haveIntersect = 1;
                        intersectPoint = [ intersectPoint ; L2B ];
                    end
                    
                    % Checking if point L1A is an intersection or not
                    arcCos_Theta_L1A_1 = dot( L2B-L2A , L1A-L2A ) / ( norm(L2B-L2A) * norm(L1A-L2A) );
                    arcCos_Theta_L1A_2 = dot( L2A-L2B , L1A-L2B ) / ( norm(L2A-L2B) * norm(L1A-L2B) );
                    if ( arcCos_Theta_L1A_1 * arcCos_Theta_L1A_2 > Epsilon ) && ( intersectNr <= 2 )
                        intersectNr = intersectNr + 1;
                        Geostatus.haveIntersect = 1;
                        intersectPoint = [ intersectPoint ; L1A ];
                    end
                    
                    % Checking if point L1B is an intersection or not
                    arcCos_Theta_L1B_1 = dot( L2B-L2A , L1B-L2A ) / ( norm(L2B-L2A) * norm(L1B-L2A) );
                    arcCos_Theta_L1B_2 = dot( L2A-L2B , L1B-L2B ) / ( norm(L2A-L2B) * norm(L1B-L2B) );
                    if ( arcCos_Theta_L1B_1 * arcCos_Theta_L1B_2 > Epsilon ) && ( intersectNr <= 2 )
                        intersectNr = intersectNr + 1;
                        Geostatus.haveIntersect = 1;
                        intersectPoint = [ intersectPoint ; L1B ];
                    end
                end
                
            else
                % The lines are either intersected or skew.
                Geostatus.areParallel  = 0;
                Geostatus.areCollinear = 0;
                t = det( [ (L2A-L1A)', V2' , cross(V1,V2)' ] ) / norm(cross(V1,V2))^2;
                s = det( [ (L2A-L1A)', V1' , cross(V1,V2)' ] ) / norm(cross(V1,V2))^2;
                
                if norm( (L1A + t*V1) - (L2A + s*V2) ) < Epsilon
                    % Intersection occurs
                    Geostatus.areSkew = 0;
                    if ( t >= Epsilon && t <= 1 + Epsilon ) && ( s >= Epsilon && s <= 1 + Epsilon )
                        % Intersection point lies in the line segments
                        Geostatus.haveIntersect  = 1;
                        intersectPoint = L1A + V1.*t;
                    else
                        % Intersection point lies out of the line segments
                        Geostatus.haveIntersect  = 0;
                    end
                else
                    % The lines are skew
                    Geostatus.haveIntersect  = 0;
                    Geostatus.areSkew = 1;
                end
            end
        end
        %
        function IsPointOnTheLineSegment = Is_Point_One_LineSegment(obj, Point, Epsilon)
            IsPointOnTheInfiniteLine = obj.Is_Point_On_InfiniteLine(Point, Epsilon);
            if IsPointOnTheInfiniteLine && dot( (obj.PointA - Point) , (obj.PointB - Point) ) < 0
                IsPointOnTheLineSegment = 1;
            else
                IsPointOnTheLineSegment = 0;
            end
        end
    end
end