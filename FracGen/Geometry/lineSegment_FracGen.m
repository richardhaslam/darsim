% class lineSegment_FracGen for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim 2 Reservoir Simulator
% Author: Mousa HosseiniMehr
% TU Delft
% Created: 2018-04-20
% Last modified: 2019-02-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef lineSegment_FracGen < lineInf_FracGen
    properties
        PointA
        PointB
        PointM
        Points
        AB_vec
    end
    methods
        function obj = lineSegment_FracGen(PointA,PointB)
            %%
            obj.InitializeLineInf(PointB-PointA);
            obj.PointA = PointA;
            obj.PointB = PointB;
            obj.AB_vec = PointB - PointA;
        end
        function AddPointM(obj)
            obj.PointM = ( obj.PointA + obj.PointB ) / 2;
            obj.Point0 = obj.PointM;
        end
        function [Geostatus, intersectPoint] = Obtain_LineSegment_LineSengment_Intersection(obj, Line, Epsilon)
            %%
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
            if ( norm( cross(obj.unit_vec, Line.unit_vec) ) < Epsilon )
                % The lines are either parallel or even collinear.
                Geostatus.areParallel = 1;
                Geostatus.areSkew = 0;
                
                % checking if there is any distance between the prarallel line segment
                parallelDistance = norm( cross( (Line.PointA - obj.PointA) , obj.unit_vec ) );
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
                        intersectPoint = [ intersectPoint , L2A ];
                    end
                    
                    % Checking if point L2B is an intersection or not
                    arcCos_Theta_L2B_1 = dot( L1B-L1A , L2B-L1A ) / ( norm(L1B-L1A) * norm(L2B-L1A) );
                    arcCos_Theta_L2B_2 = dot( L1A-L1B , L2B-L1B ) / ( norm(L1A-L1B) * norm(L2B-L1B) );
                    if ( arcCos_Theta_L2B_1 * arcCos_Theta_L2B_2 > Epsilon ) && ( intersectNr <= 2 )
                        intersectNr = intersectNr + 1;
                        Geostatus.haveIntersect = 1;
                        intersectPoint = [ intersectPoint , L2B ];
                    end
                    
                    % Checking if point L1A is an intersection or not
                    arcCos_Theta_L1A_1 = dot( L2B-L2A , L1A-L2A ) / ( norm(L2B-L2A) * norm(L1A-L2A) );
                    arcCos_Theta_L1A_2 = dot( L2A-L2B , L1A-L2B ) / ( norm(L2A-L2B) * norm(L1A-L2B) );
                    if ( arcCos_Theta_L1A_1 * arcCos_Theta_L1A_2 > Epsilon ) && ( intersectNr <= 2 )
                        intersectNr = intersectNr + 1;
                        Geostatus.haveIntersect = 1;
                        intersectPoint = [ intersectPoint , L1A ];
                    end
                    
                    % Checking if point L1B is an intersection or not
                    arcCos_Theta_L1B_1 = dot( L2B-L2A , L1B-L2A ) / ( norm(L2B-L2A) * norm(L1B-L2A) );
                    arcCos_Theta_L1B_2 = dot( L2A-L2B , L1B-L2B ) / ( norm(L2A-L2B) * norm(L1B-L2B) );
                    if ( arcCos_Theta_L1B_1 * arcCos_Theta_L1B_2 > Epsilon ) && ( intersectNr <= 2 )
                        intersectNr = intersectNr + 1;
                        Geostatus.haveIntersect = 1;
                        intersectPoint = [ intersectPoint , L1B ];
                    end
                end
                
            else
                % The lines are either intersected or skew.
                Geostatus.areParallel  = 0;
                Geostatus.areCollinear = 0;
                t = det( [ (L2A-L1A), V2 , cross(V1,V2) ] ) / norm(cross(V1,V2))^2;
                s = det( [ (L2A-L1A), V1 , cross(V1,V2) ] ) / norm(cross(V1,V2))^2;
                
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
    end
end