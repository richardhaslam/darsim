% Class of infinite plane for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2020-03-09
%Last modified: 2020-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef plane_DARSim < handle
    properties
        n_vec
        Eq
        isParallelToRefSurface
    end
    methods
        function obj = plane_DARSim(n_vec,point)
            %% Eq of the plane follows ax+by+cz=d.
            if nargin==2
                obj.n_vec = n_vec / norm(n_vec);
                obj.Eq.a = obj.n_vec(1);
                obj.Eq.b = obj.n_vec(2);
                obj.Eq.c = obj.n_vec(3);
                obj.Eq.d = obj.Eq.a * point(1) + obj.Eq.b * point(2) + obj.Eq.c * point(3);
            end
        end
        function InitializePlane(obj,n_vec,point)
            %% Eq of the plane follows ax+by+cz=d.
            obj.n_vec = n_vec / norm(n_vec);
            obj.Eq.a = obj.n_vec(1);
            obj.Eq.b = obj.n_vec(2);
            obj.Eq.c = obj.n_vec(3);
            obj.Eq.d = obj.Eq.a * point(1) + obj.Eq.b * point(2) + obj.Eq.c * point(3);
        end
        function CheckIfParallelToRefSurface(obj,Epsilon)
            %%
            if     ( abs(obj.Eq.b) < Epsilon ) && ( abs(obj.Eq.c) < Epsilon ), obj.isParallelToRefSurface = "AlongYZ";
            elseif ( abs(obj.Eq.a) < Epsilon ) && ( abs(obj.Eq.c) < Epsilon ), obj.isParallelToRefSurface = "AlongXZ";
            elseif ( abs(obj.Eq.a) < Epsilon ) && ( abs(obj.Eq.b) < Epsilon ), obj.isParallelToRefSurface = "AlongXY";
            else                                                                         , obj.isParallelToRefSurface = "None";
            end
        end
        function [Geostatus, intersectPoint] = Obtain_Plane_Line_Intersection(obj, Line, Epsilon)
            %%
            Geostatus.areParallel = NaN;
            Geostatus.areCoplanar = NaN;
            Geostatus.haveIntersect = NaN;
            intersectPoint = [];
            
            % Obtaining the geostatus between the line and the plane
            if abs( dot(obj.n_vec,Line.unit_vec) ) < Epsilon
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
                intersectPoint = [x;y;z];
            end
        end
        function [Geostatus, intersectLine] = Obtian_Plane_Plane_Intersection(obj, Plane, Epsilon)
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
                % d2 * (a1/a2) - d1 = 0 or d2 * (b1/b2) - d1 = 0 or d2 * (c1/c2) - d1 = 0
                if abs( Plane.Equation.d * (obj.equation.a/Plane.Equation.a) - obj.equation.d ) < Epsilon || ...
                   abs( Plane.Equation.d * (obj.equation.b/Plane.Equation.b) - obj.equation.d ) < Epsilon || ...
                   abs( Plane.Equation.d * (obj.equation.c/Plane.Equation.c) - obj.equation.d ) < Epsilon
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
                intersectLine = line_FracGen(unit_vec,intersectPoint);   
            end
        end
    end
end