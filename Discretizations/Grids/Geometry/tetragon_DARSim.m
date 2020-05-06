% Class of polygon for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2020-03-09
%Last modified: 2020-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef tetragon_DARSim < polygon_DARSim
    properties
        PointA
        PointB
        PointC
        PointD
        PointM
        Points
        AB_vec
        AD_vec
        AC_vec
        BC_vec
        CD_vec
        DA_vec
        Length_AB
        Length_BC
        Length_CD
        Length_DA
        AllCornersAreCoplanar
    end
    methods
        function obj = tetragon_DARSim(pointA,pointB,pointC,pointD,pos)
            %%
            obj.N_Vertex = 4;
            if nargin>=4
                obj.InitializeTetragon(pointA,pointB,pointC,pointD);
            end
            if nargin==5
                obj.AddPosition(pos);
            end
        end
        function InitializeTetragon(obj,pointA,pointB,pointC,pointD,pos)
            %%
            obj.AB_vec = pointB - pointA;
            obj.AD_vec = pointD - pointA;
            n_vec = cross(obj.AB_vec,obj.AD_vec);
            pointM = (pointA+pointB+pointC+pointD)/4;
            obj.InitializePlane(n_vec,pointM);
            obj.PointA = pointA;
            obj.PointB = pointB;
            obj.PointC = pointC;
            obj.PointD = pointD;
            obj.PointM = pointM;
            obj.Points = [pointA,pointB,pointC,pointD,pointM];
            obj.AC_vec = obj.PointC - obj.PointA;
            obj.BC_vec = obj.PointC - obj.PointB;
            obj.CD_vec = obj.PointD - obj.PointC;
            obj.DA_vec = obj.PointA - obj.PointD;
            obj.Check_If_All_Corners_Coplanar();
            if nargin==6
                obj.AddPosition(pos);
            end
        end
        function Check_If_All_Corners_Coplanar(obj)
            % From pointA, we draw vectors to each of the points B,C and D.
            % Assume the vectors are v1, v2 and v3. The four points are
            % coplanar if the volume created by these three vectors is zero.
            % Namely: v1 . (v2 x v3) = 0;
            dx = norm(obj.PointA - obj.PointC);
            dy = norm(obj.PointB - obj.PointD);
            Epsilon = 1e-10 * mean([dx,dy]);
            if abs( dot( obj.AB_vec, cross(obj.AC_vec, obj.AD_vec) ) ) > Epsilon
                % warning('This is not a tetragon, because the four corners are not coplanar!');
                obj.AllCornersAreCoplanar = 0;
            else
                obj.AllCornersAreCoplanar = 1;
            end
        end
        function isInside = Is_Point_Inside_Tetrahedron(obj, point , Epsilon)
            %%
            isInside  = NaN;
            if dot( obj.n_vec , (point-obj.PointM) ) > Epsilon
                % The Point is not on the plane
                isInside = 0;
                return;
            end
            
            Line1 = lineSegment_FracGen(obj.PointM , point);
            
            Line2 = lineSegment_FracGen(obj.PointA , obj.PointB);
            [Geostatus, IntersectPoint] = Line1.Obtain_LineSegment_LineSengment_Intersection( Line2, Epsilon );
            if     ( Geostatus.haveIntersect == 0 )             ,  isInside = 1;
            elseif ( norm (IntersectPoint - point ) < Epsilon ) ,  isInside = 1;
            else                                                ,  isInside = 0;  return;  end
            
            Line2 = lineSegment_FracGen(obj.PointB , obj.PointC);
            [Geostatus, IntersectPoint] = Line1.Obtain_LineSegment_LineSengment_Intersection( Line2, Epsilon );
            if     ( Geostatus.haveIntersect == 0 )             ,  isInside = 1;
            elseif ( norm (IntersectPoint - point ) < Epsilon ) ,  isInside = 1;
            else                                                ,  isInside = 0;  return;  end
            
            Line2 = lineSegment_FracGen(obj.PointC , obj.PointD);
            [Geostatus, IntersectPoint] = Line1.Obtain_LineSegment_LineSengment_Intersection( Line2, Epsilon );
            if     ( Geostatus.haveIntersect == 0 )             ,  isInside = 1;
            elseif ( norm (IntersectPoint - point ) < Epsilon ) ,  isInside = 1;
            else                                                ,  isInside = 0;  return;  end
            
            Line2 = lineSegment_FracGen(obj.PointD , obj.PointA);
            [Geostatus, IntersectPoint] = Line1.Obtain_LineSegment_LineSengment_Intersection( Line2, Epsilon );
            if     ( Geostatus.haveIntersect == 0 )             ,  isInside = 1;
            elseif ( norm (IntersectPoint - point ) < Epsilon ) ,  isInside = 1;
            else                                                ,  isInside = 0;  return;  end
        end
        function [Geostatus, IntersectPoints] = Obtain_Tetrahedron_Tetrahedron_Intersection(obj, tetrahedron, Epsilon )
            %%
            Geostatus.areParallel   = NaN;
            Geostatus.areCoplanar   = NaN;
            Geostatus.haveIntersect = NaN;
            IntersectPoints         = [];
            
            % Obtaining the geostatus between the planes
            if norm( cross( obj.n_vec , tetrahedron.n_vec ) ) < Epsilon
                % The planes are parallel
                Geostatus.areParallel = 1;
                % If the equations of both planes are multiple of each
                % other, then they are coplanar:
                % To check this, we put a point from the 1st plane into
                % 2nd plane's equation. If the equation holds, it
                % means that the point lies within both planes resulting in
                % coplanarity of these two planes.
                % if a2*x1 + b2*y1 + c2*z1 = d2, then the planes are coplanar.
                if abs( Plane.Equation.a * obj.PointM(1) + Plane.Equation.b * obj.PointM(2) + ...
                        Plane.Equation.c * obj.PointM(3) - Plane.Equation.d ) < Epsilon
                    % The planes are coplanar
                    Geostatus.areCoplanar = 1;
                    
                    % Checking the intersections between each two edges of plane segments
                    for i = [1,2,3,4]
                        intersectNr = 0;
                        for j = [2,3,4,1]
                            Line1 = lineSegment_DARSim(obj.Points(:,i)  ,obj.Points(:,j)  );
                            Line2 = lineSegment_DARSim(tetrahedron.Points(:,i),tetrahedron.Points(:,j));
                            [lineGeostatus, lineIntersectPoint] = Line1.Obtain_LineSegment_LineSengment_Intersection(Line2, Epsilon);
                            if lineGeostatus.haveIntersect == 1
                                Geostatus.haveIntersect = 1;
                                intersectNr = intersectNr + 1;
                                IntersectPoints = [IntersectPoints , lineIntersectPoint];
                            end
                            if intersectNr == 2,  continue;  end  % each edge of 1st plane segment can have intersection with maximum two edges of the 2nd plane segment
                        end
                    end
                    
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
                
                % Obtaining the unit vector of the intersection line
                intL_Vec = cross( obj.n_vec , tetrahedron.n_vec );
                intL_Vec = intL_Vec / norm(intL_Vec);
                
                % Obtaining a Point on the intersection line with Z=0 if possible
                P2M = (tetrahedron.PointA + tetrahedron.PointC )/2;
                Axy = [ obj.n_vec(1) , obj.n_vec(2) ; tetrahedron.n_vec(1)   , tetrahedron.n_vec(2) ];
                RHS = [ obj.Equation.d   - obj.n_vec(3)*P2M(3)
                    tetrahedron.Equation.d - tetrahedron.n_vec(3)*P2M(3) ];
                if abs(det(Axy)) > Epsilon
                    intL_Point0    = zeros(3,1);
                    Unknowns       = Axy \ RHS;
                    intL_Point0(1) = Unknowns(1);
                    intL_Point0(2) = Unknowns(2);
                    intL_Point0(3) = P2M(3);
                else
                    % Or, obtaining a Point on the intersection line with Y=0 if possible
                    Axz = [ obj.n_vec(1) , obj.n_vec(3) ; tetrahedron.n_vec(1)  , tetrahedron.n_vec(3) ];
                    RHS = [ obj.Equation.d   - obj.n_vec(2)*P2M(2)
                        tetrahedron.Equation.d - tetrahedron.n_vec(2)*P2M(2) ];
                    if abs(det(Axz)) > Epsilon
                        intL_Point0    = zeros(3,1);
                        Unknowns       = Axz \ RHS;
                        intL_Point0(1) = Unknowns(1);
                        intL_Point0(3) = Unknowns(2);
                        intL_Point0(2) = P2M(2);
                    else
                        % Or, obtaining a Point on the intersection line with X=0 if possible
                        Ayz = [ obj.n_vec(2) , obj.n_vec(3) ; tetrahedron.n_vec(2)  , tetrahedron.n_vec(3) ];
                        RHS = [ obj.Equation.d   - obj.n_vec(1)*P2M(1)
                            tetrahedron.Equation.d - tetrahedron.n_vec(1)*P2M(1) ];
                        intL_Point0    = zeros(3,1);
                        Unknowns       = Ayz \ RHS;
                        intL_Point0(2) = Unknowns(1);
                        intL_Point0(3) = Unknowns(2);
                        intL_Point0(1) = P2M(1);
                    end
                end
                
                % Assuming two end-Points for the intersection line
                intL_A = intL_Point0 - intL_Vec * max( norm(obj.AC_vec) , norm(tetrahedron.AC_vec) ) * 1e3;
                intL_B = intL_Point0 + intL_Vec * max( norm(obj.AC_vec) , norm(tetrahedron.AC_vec) ) * 1e3;
                intersectionLine = lineSegment_FracGen(intL_A  ,intL_B);
                
                % Checking the intersection between the obtained line segment and each side of each plane segment
                allPoints = [obj.Points(:,1:4), tetrahedron.Points(:,1:4)];
                ind1 = [1,2,3,4,5,6,7,8];
                ind2 = [2,3,4,1,6,7,8,5];
                intersectNr = 0;
                for i = 1 : 8
                    Line_temp = lineSegment_FracGen( allPoints(:,ind1(i)) , allPoints(:,ind2(i)) );
                    [lineGeostatus, lineIntersectPoint] = Line_temp.Obtain_LineSegment_LineSengment_Intersection( intersectionLine, Epsilon );
                    if lineGeostatus.haveIntersect == 1
                        Geostatus.haveIntersect = 1;
                        intersectNr = intersectNr + 1;
                        IntersectPoints = [IntersectPoints , lineIntersectPoint];
                    end
                    if intersectNr == 4,  continue;  end % intersection line can have intersection with maximum two edges of each plane segment (four intersections max)
                end
                
            end
            
        end
        function AvgDist = Obtain_Average_Dsitance_Line_From_Tetrahedron(obj, Line, Refinement)
            % This functions calculates the average distance between a
            % tetrahedron and a line. The line should lie within the
            % plane surface.
            
            % PlaneEquation has 4 set of values of "a,b,c,d" and follows the equation
            % pattern "ax+by+cz=d".
            
            % "Refinement" is number of devisions in every direction. For example,
            % if the refinement number is 7, the plane segment is devided into 7^2 = 49 smaller
            % plane segments to calculate the average distance.
            
            Refinement = round(max(Refinement, 2));
            AvgDist = 0;
 
            AD_Coords = [linsapce(obj.PointA(1), obj.PointD(1), Refinement+1)
                         linsapce(obj.PointA(2), obj.PointD(2), Refinement+1)
                         linsapce(obj.PointA(3), obj.PointD(3), Refinement+1)];
                     
            BC_Coords = [linsapce(obj.PointB(1), obj.PointC(1), Refinement+1)
                         linsapce(obj.PointB(2), obj.PointC(2), Refinement+1)
                         linsapce(obj.PointB(3), obj.PointC(3), Refinement+1)];
            
            GridCoords  = zeros( Refinement+1, Refinement+1, 3 );
            for i = 1 : Refinement+1
                GridCoords(:,i,1) = linspace( AD_Coords(i,1) , BC_Coords(i,1) , Refinement+1 );                          
                GridCoords(:,i,2) = linspace( AD_Coords(i,2) , BC_Coords(i,2) , Refinement+1 );
                GridCoords(:,i,3) = linspace( AD_Coords(i,3) , BC_Coords(i,3) , Refinement+1 );
            end
            
            x1 = Line.Point0  - Line.unit_vec * 1e3;
            x2 = Line.Point0  + Line.unit_vec * 1e3;
            for j = 1 : Refinement+1
                for i = 1 : Refinement+1
                    x0 = [GridCoords(i,j,1);GridCoords(i,j,2);GridCoords(i,j,3)];
                    AvgDist = AvgDist + norm(cross(x2-x1,x1-x0))/norm(x2-x1);
                end
            end
            AvgDist = AvgDist / (Refinement+1)^2 ;
        end
    end
end