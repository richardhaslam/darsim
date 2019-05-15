% Class lineInf_FracGen for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim 2 Reservoir Simulator
% Author: Mousa HosseiniMehr
% TU Delft
% Created: 2018-04-20
% Last modified: 2019-02-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef lineInf_FracGen < handle
    properties
        Point0
        unit_vec
        Equation % The equation of line is applied as follows x = a*t+x0; y = b*t+y0; z = c*t+z0;
        isParallelToRefAxis
        isParallelToRefSurface
    end
    methods
        function obj = lineInf_FracGen(unit_vec)
            %%
            if nargin==1
                obj.unit_vec = unit_vec / norm(unit_vec);
            end
        end
        function InitializeLineInf(obj,unit_vec)
            %%
            obj.unit_vec = unit_vec / norm(unit_vec);
        end
        function AddPoint0(obj,Point0)
            obj.Point0 = Point0;
        end
        function AddEquation(obj,Point0)
            %
            if nargin==2 && isempty(obj.Point0)
                obj.AddPoint0(Point0);
            end
            obj.Equation.a = obj.unit_vec(1);
            obj.Equation.b = obj.unit_vec(2);
            obj.Equation.c = obj.unit_vec(3);
            obj.Equation.x0 = obj.Point0(1);
            obj.Equation.y0 = obj.Point0(2);
            obj.Equation.z0 = obj.Point0(3);
        end
        function [Geostatus, intersectPoint] = Obtain_LineInf_LineInf_Intersection(obj, Line, Epsilon)
            %%
            Geostatus.areParallel = NaN;
            Geostatus.areCollinear = NaN;
            Geostatus.haveIntersect = NaN;
            intersectPoint = [];
            
            % Obtaining the geostatus between the lines
            if ( norm( cross(obj.unit_vec, Line.unit_vec) ) < Epsilon )
                % The lines are either parallel or even collinear.
                Geostatus.areParallel = 1;
                
                % checking if there is any distance between the prarallel line segment
                parallelDistance = norm( cross( (Line.Point0 - obj.Point0) , obj.unit_vec ) );
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
                L1A = obj.Point0  - obj.unit_vec * 1e3;
                L1B = obj.Point0  + obj.unit_vec * 1e3;
                L2A = Line.Point0 - Line.unit_vec * 1e3;
                L2B = Line.Point0 + Line.unit_vec * 1e3;
                V1 = L1B - L1A;
                V2 = L2B - L2A;
                t = det( [ (L2A-L1A), V2 , cross(V1,V2) ] ) / norm(cross(V1,V2))^2;
                s = det( [ (L2A-L1A), V1 , cross(V1,V2) ] ) / norm(cross(V1,V2))^2;
                
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
    end
end