% class plate for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2018-04-20
%Last modified: 2018-04-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef planeSegment_FracGen < planeInf_FracGen
    properties
        PointA
        PointB
        PointC
        PointD
        PointM
        Points
        Length_AB
        Width_AD
        AB_vec
        AD_vec
        AC_vec
        BC_vec
        CD_vec
        DA_vec
    end
    methods
        function obj = planeSegment_FracGen(PointA,PointB,PointC,PointD)
            %%
            if nargin==4
                AB_vec = PointB - PointA;
                AD_vec = PointD - PointA;
                n_vec = cross( AB_vec , AD_vec );
                PointM = (PointB + PointD ) /2;
                obj.InitializePlaneInf(n_vec,PointM);
                obj.PointA = PointA;
                obj.PointB = PointB;
                obj.PointC = PointC;
                obj.PointD = PointD;
                obj.PointM = PointM;
                obj.Points = [PointA,PointB,PointC,PointD,PointM];
                obj.AB_vec = AB_vec;
                obj.AD_vec = AD_vec;
                obj.AC_vec = obj.PointC - obj.PointA;
                obj.BC_vec = obj.PointC - obj.PointB;
                obj.CD_vec = obj.PointD - obj.PointC;
                obj.DA_vec = obj.PointA - obj.PointD;
            end
        end
        function InitializePlaneSegment(obj,PointA,PointB,PointC,PointD)
            %%
            AB_vec = PointB - PointA;
            AD_vec = PointD - PointA;
            n_vec = cross( AB_vec , AD_vec );
            PointM = (PointB + PointD ) /2;
            obj.InitializePlaneInf(n_vec,PointM);
            obj.PointA = PointA;
            obj.PointB = PointB;
            obj.PointC = PointC;
            obj.PointD = PointD;
            obj.PointM = PointM;
            obj.Points = [PointA,PointB,PointC,PointD,PointM];
            obj.AB_vec = AB_vec;
            obj.AD_vec = AD_vec;
            obj.AC_vec = obj.PointC - obj.PointA;
            obj.BC_vec = obj.PointC - obj.PointB;
            obj.CD_vec = obj.PointD - obj.PointC;
            obj.DA_vec = obj.PointA - obj.PointD;
        end
        function isInside = Is_Point_Inside_PlaneSegment( obj, Point , Epsilon)
            %%
            isInside  = NaN;
            if dot( obj.n_vec , (Point-obj.PointM) ) > Epsilon
                % The Point is not on the plane
                isInside = 0;
                return;
            end
            
            Line1 = lineSegment_FracGen(obj.PointM , Point);
            
            Line2 = lineSegment_FracGen(obj.PointA , obj.PointB);
            [Geostatus, IntersectPoint] = Line1.Obtain_LineSegment_LineSengment_Intersection( Line2, Epsilon );
            if     ( Geostatus.haveIntersect == 0 )             ,  isInside = 1;
            elseif ( norm (IntersectPoint - Point ) < Epsilon ) ,  isInside = 1;
            else                                                ,  isInside = 0;  return;  end
            
            Line2 = lineSegment_FracGen(obj.PointB , obj.PointC);
            [Geostatus, IntersectPoint] = Line1.Obtain_LineSegment_LineSengment_Intersection( Line2, Epsilon );
            if     ( Geostatus.haveIntersect == 0 )             ,  isInside = 1;
            elseif ( norm (IntersectPoint - Point ) < Epsilon ) ,  isInside = 1;
            else                                                ,  isInside = 0;  return;  end
            
            Line2 = lineSegment_FracGen(obj.PointC , obj.PointD);
            [Geostatus, IntersectPoint] = Line1.Obtain_LineSegment_LineSengment_Intersection( Line2, Epsilon );
            if     ( Geostatus.haveIntersect == 0 )             ,  isInside = 1;
            elseif ( norm (IntersectPoint - Point ) < Epsilon ) ,  isInside = 1;
            else                                                ,  isInside = 0;  return;  end
            
            Line2 = lineSegment_FracGen(obj.PointD , obj.PointA);
            [Geostatus, IntersectPoint] = Line1.Obtain_LineSegment_LineSengment_Intersection( Line2, Epsilon );
            if     ( Geostatus.haveIntersect == 0 )             ,  isInside = 1;
            elseif ( norm (IntersectPoint - Point ) < Epsilon ) ,  isInside = 1;
            else                                                ,  isInside = 0;  return;  end
        end
        function [Geostatus, IntersectPoints] = Obtain_PlaneSegment_PlaneSegment_Intersection(obj, Plane, Epsilon )
            %%
            Geostatus.areParallel   = NaN;
            Geostatus.areCoplanar   = NaN;
            Geostatus.haveIntersect = NaN;
            IntersectPoints         = [];
            
            % Obtaining the geostatus between the planes
            if norm( cross( obj.n_vec , Plane.n_vec ) ) < Epsilon
                % The planes are parallel
                Geostatus.areParallel = 1;
                % If the equations of both planes are multiple of each
                % other, then they are coplanar:
                % d2 * (a1/a2) - d1 = 0
                if abs( Plane.Equation.d * (obj.Equation.a/Plane.Equation.a) - obj.Equation.d ) < Epsilon
                    % The planes are coplanar
                    Geostatus.areCoplanar = 1;
                    
                    % Checking the intersections between each two edges of plane segments
                    for i = [1,2,3,4]
                        intersectNr = 0;
                        for j = [2,3,4,1]
                            Line1 = lineSegment_FracGen(obj.Points(:,i)  ,obj.Points(:,j)  );
                            Line2 = lineSegment_FracGen(Plane.Points(:,i),Plane.Points(:,j));
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
                intL_Vec = cross( obj.n_vec , Plane.n_vec );
                intL_Vec = intL_Vec / norm(intL_Vec);
                
                % Obtaining a Point on the intersection line with Z=0 if possible
                P2M = (Plane.PointA + Plane.PointC )/2;
                Axy = [ obj.n_vec(1) , obj.n_vec(2) ; Plane.n_vec(1)   , Plane.n_vec(2) ];
                RHS = [ obj.Equation.d   - obj.n_vec(3)*P2M(3)
                    Plane.Equation.d - Plane.n_vec(3)*P2M(3) ];
                if abs(det(Axy)) > Epsilon
                    intL_Point0    = zeros(3,1);
                    Unknowns       = Axy \ RHS;
                    intL_Point0(1) = Unknowns(1);
                    intL_Point0(2) = Unknowns(2);
                    intL_Point0(3) = P2M(3);
                else
                    % Or, obtaining a Point on the intersection line with Y=0 if possible
                    Axz = [ obj.n_vec(1) , obj.n_vec(3) ; Plane.n_vec(1)  , Plane.n_vec(3) ];
                    RHS = [ obj.Equation.d   - obj.n_vec(2)*P2M(2)
                        Plane.Equation.d - Plane.n_vec(2)*P2M(2) ];
                    if abs(det(Axz)) > Epsilon
                        intL_Point0    = zeros(3,1);
                        Unknowns       = Axz \ RHS;
                        intL_Point0(1) = Unknowns(1);
                        intL_Point0(3) = Unknowns(2);
                        intL_Point0(2) = P2M(2);
                    else
                        % Or, obtaining a Point on the intersection line with X=0 if possible
                        Ayz = [ obj.n_vec(2) , obj.n_vec(3) ; Plane.n_vec(2)  , Plane.n_vec(3) ];
                        RHS = [ obj.Equation.d   - obj.n_vec(1)*P2M(1)
                            Plane.Equation.d - Plane.n_vec(1)*P2M(1) ];
                        intL_Point0    = zeros(3,1);
                        Unknowns       = Ayz \ RHS;
                        intL_Point0(2) = Unknowns(1);
                        intL_Point0(3) = Unknowns(2);
                        intL_Point0(1) = P2M(1);
                    end
                end
                
                % Assuming two end-Points for the intersection line
                intL_A = intL_Point0 - intL_Vec * max( norm(obj.AC_vec) , norm(Plane.AC_vec) ) * 1e3;
                intL_B = intL_Point0 + intL_Vec * max( norm(obj.AC_vec) , norm(Plane.AC_vec) ) * 1e3;
                intersectionLine = lineSegment_FracGen(intL_A  ,intL_B);
                
                % Checking the intersection between the obtained line segment and each side of each plane segment
                allPoints = [obj.Points(:,1:4), Plane.Points(:,1:4)];
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
        function AvgDist = Obtain_Average_Dsitance_LineInf_From_PlaneSegment(obj, Line, Refinement)
            % This functions calculates the average distance between a
            % plane segment and a line. The line should lie within the
            % plane surface.
            
            % PlaneEquation has 4 set of values of "a,b,c,d" and follows the equation
            % pattern "ax+by+cz=d".
            
            % "Refinement" is number of devisions in every direction. For example,
            % if the refinement number is 7, the plane segment is devided into 7^2 = 49 smaller
            % plane segments to calculate the average distance.
            
            Refinement = round(max(Refinement, 2));
            
            AvgDist = 0;
            
            obj.Length_AB = norm( obj.AB_vec );
            obj.Width_AD  = norm( obj.AD_vec );
            D_Length_AB_vec = obj.AB_vec / Refinement;
            % D_Width_AD_vec  = obj.AD_vec / Refinement;
            
            GridCoords  = zeros( Refinement+1, Refinement+1, 3 );
            for i = 1 : Refinement+1
                GridCoords(i,:,1) = linspace( obj.PointA(1) + (i-1)*D_Length_AB_vec(1) , obj.PointD(1) + (i-1)*D_Length_AB_vec(1) , Refinement+1 );                          
                GridCoords(i,:,2) = linspace( obj.PointA(2) + (i-1)*D_Length_AB_vec(2) , obj.PointD(2) + (i-1)*D_Length_AB_vec(2) , Refinement+1 );                        
                GridCoords(i,:,3) = linspace( obj.PointA(3) + (i-1)*D_Length_AB_vec(3) , obj.PointD(3) + (i-1)*D_Length_AB_vec(3) , Refinement+1 );
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