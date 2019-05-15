% class plate for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2018-04-20
%Last modified: 2018-04-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef planeInf_FracGen < handle
    properties
        n_vec
        Equation
        isParallelToRefSurface
    end
    methods
        function obj = planeInf_FracGen(n_vec,Point)
            %% Equation of the plane follows ax+by+cz=d.
            if nargin==2
                obj.n_vec = n_vec / norm(n_vec);
                obj.Equation.a = obj.n_vec(1);
                obj.Equation.b = obj.n_vec(2);
                obj.Equation.c = obj.n_vec(3);
                obj.Equation.d = obj.Equation.a * Point(1) + obj.Equation.b * Point(2) + obj.Equation.c * Point(3);
            end
        end
        function InitializePlaneInf(obj,n_vec,Point)
            %%
            obj.n_vec = n_vec / norm(n_vec);
            obj.Equation.a = obj.n_vec(1);
            obj.Equation.b = obj.n_vec(2);
            obj.Equation.c = obj.n_vec(3);
            obj.Equation.d = obj.Equation.a * Point(1) + obj.Equation.b * Point(2) + obj.Equation.c * Point(3);
        end
        function CheckIfParallelToRefSurface(obj,Epsilon)
            %%
            if     ( abs(obj.Equation.b) < Epsilon ) && ( abs(obj.Equation.c) < Epsilon ), obj.isParallelToRefSurface = "AlongYZ";
            elseif ( abs(obj.Equation.a) < Epsilon ) && ( abs(obj.Equation.c) < Epsilon ), obj.isParallelToRefSurface = "AlongXZ";
            elseif ( abs(obj.Equation.a) < Epsilon ) && ( abs(obj.Equation.b) < Epsilon ), obj.isParallelToRefSurface = "AlongXY";
            else                                                                         , obj.isParallelToRefSurface = "None";
            end
        end
        function [Geostatus, intersectPoint] = Obtain_LineInf_PlaneInf_Intersection(obj, Line, Epsilon)
            %%
            Geostatus.areParallel = NaN;
            Geostatus.areCoplanar = NaN;
            Geostatus.haveIntersect = NaN;
            intersectPoint = [];
            
            % Obtaining the geostatus between the line and the plane
            if abs( dot(obj.n_vec,Line.unit_vec) ) < Epsilon
                % The line and the plane are parallel
                Geostatus.areParallel = 1;
                if abs( obj.Equation.a*Line.Equation.x0 + obj.Equation.b*Line.Equation.y0 + obj.Equation.c*Line.Equation.z0 - obj.Equation.d ) < Epsilon
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
                
                % Line equation is: x=a't+x0, y=b't+y0, z=c't+z0.
                % Plane equation is: ax+by+cz=d.
                % By substituting line equation into plane equation, the
                % parameter "t" is obtained:
                % a(a't+x0)+b(b't+y0)+c(c't+z0)=d
                % ===> t = ( d-ax0-by0-cz0 ) / ( a.a' + b.b' + c.c' );
                % Afterwards, by putting the value of "t" parameter in the
                % line equation, the [x,y,z] of the intersection point is
                % obtained.
                
                t = ( obj.Equation.d - obj.Equation.a*Line.Equation.x0 - obj.Equation.b*Line.Equation.y0 - obj.Equation.c*Line.Equation.z0 ) ...
                    / ( obj.Equation.a*Line.Equation.a' + obj.Equation.b*Line.Equation.b' + obj.Equation.c*Line.Equation.c' );
                
                x = Line.Equation.a * t + Line.Equation.x0;
                y = Line.Equation.b * t + Line.Equation.y0;
                z = Line.Equation.c * t + Line.Equation.z0;
                intersectPoint = [x;y;z];
            end
        end
        function [Geostatus, intersectLine] = Obtian_PlaneInf_PlaneInf_Intersection(obj, Plane, Epsilon)
            %%
            Geostatus.areParallel = NaN;
            Geostatus.areCoplanar = NaN;
            Geostatus.haveIntersect = NaN;
            intersectLine = [];
            % Obtaining the geostatus between the planes
            if norm( cross( obj.n_vec , Plane.n_vec ) ) < Epsilon
                % The planes are parallel
                Geostatus.areParallel = 1;
                % If the equations of both planes are multiple of each
                % other, then they are coplanar:
                % d2 * (a1/a2) - d1 = 0
                if abs( Plane.Equation.d * (obj.equation.a/Plane.Equation.a) - obj.equation.d ) < Epsilon
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
                unit_vec = cross(obj.n_vec,Plane.n_vec);
                unit_vec = unit_vec / norm(unit_vec);
                
                % Obtaining a Point on the intersection line:
                Axy = [ obj.Equation.a , obj.Equation.b ; Plane.Equation.a , Plane.Equation.b ];
                Axz = [ obj.Equation.a , obj.Equation.c ; Plane.Equation.a , Plane.Equation.c ];
                Ayz = [ obj.Equation.b , obj.Equation.c ; Plane.Equation.b , Plane.Equation.c ];
                RHS = [ obj.Equation.d ; Plane.Equation.d ];
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
                intersectLine = line_FracGen(unit_vec,intersectPoint);   
            end
        end
    end
end