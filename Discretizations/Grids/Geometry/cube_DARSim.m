% Class of cube for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2020-03-09
%Last modified: 2020-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef cube_DARSim < hexahedron_DARSim
    properties
        PointA
        PointB
    end
    methods
        %%
        function obj = cube_DARSim(pointA,pointB)
            obj.NumOfVertex = 8;
            obj.NumOfFace = 6;
            obj.Face = cube_DARSim.empty;
            if nargin >= 2
                obj.InitializeCube(pointA,pointB);
            end
        end
        %%
        function InitializeCube(obj,pointA,pointB)
            obj.PointA = pointA;
            obj.PointB = pointB;
            NW_Top = [ pointA(1) ; pointB(2) ; pointB(3) ];
            SW_Top = [ pointA(1) ; pointA(2) ; pointB(3) ];
            SE_Top = [ pointB(1) ; pointA(2) ; pointB(3) ];
            NE_Top = [ pointB(1) ; pointB(2) ; pointB(3) ];
            NW_Bot = [ pointA(1) ; pointB(2) ; pointA(3) ];
            SW_Bot = [ pointA(1) ; pointA(2) ; pointA(3) ];
            SE_Bot = [ pointB(1) ; pointA(2) ; pointA(3) ];
            NE_Bot = [ pointB(1) ; pointB(2) ; pointA(3) ];
            obj.InitializeHexahedron(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot);
        end
        %%
        function AvgDistance = Obtain_Average_Distance_Plane_From_Cube(obj, Plane, Refinement)
            % This functions calculates the average distance between a cube and a
            % plane.
            
            % PlaneEquation has 4 set of values of "a,b,c,d" and follows the equation
            % pattern "ax+by+cz=d".
            
            % CubeCornerStart and CubeCornerEnd are the coordinates of starting and
            % ending corners of the cube (x,y,z).
            
            % "Refinement" is number of devisions in every direction. For example,
            % if the refinement number is 7, the cube is devided into 7^3 = 343 smaller
            % cubes to calculate the average distance.
            
            Refinement = round(max(Refinement, 2));
            
            AvgDistance = 0;
            
            XimCube = linspace( obj.PointA(1), obj.PointB(1), Refinement+1 );
            YimCube = linspace( obj.PointA(2), obj.PointB(2), Refinement+1 );
            ZimCube = linspace( obj.PointA(3), obj.PointB(3), Refinement+1 );
            
            a = Plane.Eq.a;
            b = Plane.Eq.b;
            c = Plane.Eq.c;
            d = Plane.Eq.d;
            for k = 1 : Refinement+1
                for j = 1 : Refinement+1
                    AvgDistance = AvgDistance + sum ( abs( a * XimCube + b * YimCube(j) + c * ZimCube(k) - d ) );
                end
            end
            
            AvgDistance = AvgDistance / sqrt( a^2 + b^2 + c^2 ) / (Refinement+1)^3 ;
        end
        %%
        function [intersectCoordFinal, intersectCoordTemp , areCoplanar] = Obtain_Cube_Tetragon_Intersection(obj, Tetragon, Epsilon)
            % Assigning reservoir properties
            % Coordinates of reservoir grid nodes in each direction
            
            XimCube = linspace( obj.PointA(1), obj.PointB(1), 2 );
            YimCube = linspace( obj.PointA(2), obj.PointB(2), 2 );
            ZimCube = linspace( obj.PointA(3), obj.PointB(3), 2 );
            
            % Initializing some variables
            doNotContinue       = 0;
            areCoplanar         = 0;
            intersectCoordTemp  = []; 
            intersectCoordFinal = [];
            
            % Matrix Cube Face #X1
            if doNotContinue == 0
                PointA = [ XimCube(1) ; YimCube(1) ; ZimCube(1) ];
                PointB = [ XimCube(1) ; YimCube(2) ; ZimCube(1) ];
                PointC = [ XimCube(1) ; YimCube(2) ; ZimCube(2) ];
                PointD = [ XimCube(1) ; YimCube(1) ; ZimCube(2) ];
                Plane_X1 = tetragon_DARSim(PointA,PointB,PointC,PointD);
                
                [Geostatus_X1, intersectPoints_X1] = Tetragon.Obtain_Tetragon_Tetragon_Intersection( Plane_X1 , Epsilon );
                if ( Geostatus_X1.areCoplanar == 1 )
                    doNotContinue      = 1;
                    areCoplanar        = 1;
                    intersectPoints_X1 = [ intersectPoints_X1 , Plane_X1.PointA , Plane_X1.PointB , Plane_X1.PointC , Plane_X1.PointD ];
                end
                intersectCoordTemp = [ intersectCoordTemp , intersectPoints_X1 ];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cube Face #X2
            if doNotContinue == 0
                PointA = [ XimCube(2) ; YimCube(1) ; ZimCube(1) ];
                PointB = [ XimCube(2) ; YimCube(2) ; ZimCube(1) ];
                PointC = [ XimCube(2) ; YimCube(2) ; ZimCube(2) ];
                PointD = [ XimCube(2) ; YimCube(1) ; ZimCube(2) ];
                Plane_X2 = tetragon_DARSim(PointA,PointB,PointC,PointD);
                
                [Geostatus_X2, intersectPoints_X2] = Tetragon.Obtain_Tetragon_Tetragon_Intersection( Plane_X2 , Epsilon );
                if ( Geostatus_X2.areCoplanar == 1 )
                    doNotContinue      = 1;
                    areCoplanar        = 1;
                    intersectPoints_X2 = [ intersectPoints_X2 , Plane_X2.PointA , Plane_X2.PointB , Plane_X2.PointC , Plane_X2.PointD ];
                end
                intersectCoordTemp = [ intersectCoordTemp , intersectPoints_X2 ];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cube Face #Y1
            if doNotContinue == 0
                PointA = [ XimCube(1) ; YimCube(1) ; ZimCube(1) ];
                PointB = [ XimCube(2) ; YimCube(1) ; ZimCube(1) ];
                PointC = [ XimCube(2) ; YimCube(1) ; ZimCube(2) ];
                PointD = [ XimCube(1) ; YimCube(1) ; ZimCube(2) ];
                Plane_Y1 = tetragon_DARSim(PointA,PointB,PointC,PointD);
                
                [Geostatus_Y1, intersectPoints_Y1] = Tetragon.Obtain_Tetragon_Tetragon_Intersection( Plane_Y1 , Epsilon );
                if ( Geostatus_Y1.areCoplanar == 1 )
                    doNotContinue      = 1;
                    areCoplanar        = 1;
                    intersectPoints_Y1 = [ intersectPoints_Y1 , Plane_Y1.PointA , Plane_Y1.PointB , Plane_Y1.PointC , Plane_Y1.PointD ];
                end
                intersectCoordTemp = [ intersectCoordTemp , intersectPoints_Y1 ];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cube Face #Y2
            if doNotContinue == 0
                PointA = [ XimCube(1) ; YimCube(2) ; ZimCube(1) ];
                PointB = [ XimCube(2) ; YimCube(2) ; ZimCube(1) ];
                PointC = [ XimCube(2) ; YimCube(2) ; ZimCube(2) ];
                PointD = [ XimCube(1) ; YimCube(2) ; ZimCube(2) ];
                Plane_Y2 = tetragon_DARSim(PointA,PointB,PointC,PointD);
                
                [Geostatus_Y2, intersectPoints_Y2] = Tetragon.Obtain_Tetragon_Tetragon_Intersection( Plane_Y2 , Epsilon );
                if ( Geostatus_Y2.areCoplanar == 1 )
                    doNotContinue      = 1;
                    areCoplanar        = 1;
                    intersectPoints_Y2 = [ intersectPoints_Y2 , Plane_Y2.PointA , Plane_Y2.PointB , Plane_Y2.PointC , Plane_Y2.PointD ];
                end
                intersectCoordTemp = [ intersectCoordTemp , intersectPoints_Y2 ];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cube Face #Z1
            if doNotContinue == 0
                PointA = [ XimCube(1) ; YimCube(1) ; ZimCube(1) ];
                PointB = [ XimCube(2) ; YimCube(1) ; ZimCube(1) ];
                PointC = [ XimCube(2) ; YimCube(2) ; ZimCube(1) ];
                PointD = [ XimCube(1) ; YimCube(2) ; ZimCube(1) ];
                Plane_Z1 = tetragon_DARSim(PointA,PointB,PointC,PointD);
                
                [Geostatus_Z1, intersectPoints_Z1] = Tetragon.Obtain_Tetragon_Tetragon_Intersection( Plane_Z1 , Epsilon );
                if ( Geostatus_Z1.areCoplanar == 1 )
                    doNotContinue      = 1;
                    areCoplanar        = 1;
                    intersectPoints_Z1 = [ intersectPoints_Z1 , Plane_Z1.PointA , Plane_Z1.PointB , Plane_Z1.PointC , Plane_Z1.PointD ];
                end
                intersectCoordTemp = [ intersectCoordTemp , intersectPoints_Z1 ];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cube Face #Z2
            if doNotContinue == 0
                PointA = [ XimCube(1) ; YimCube(1) ; ZimCube(2) ];
                PointB = [ XimCube(2) ; YimCube(1) ; ZimCube(2) ];
                PointC = [ XimCube(2) ; YimCube(2) ; ZimCube(2) ];
                PointD = [ XimCube(1) ; YimCube(2) ; ZimCube(2) ];
                Plane_Z2 = tetragon_DARSim(PointA,PointB,PointC,PointD);
                
                [Geostatus_Z2, intersectPoints_Z2] = Tetragon.Obtain_Tetragon_Tetragon_Intersection( Plane_Z2 , Epsilon );
                if ( Geostatus_Z2.areCoplanar == 1 )
                    doNotContinue = 1;
                    areCoplanar   = 1;
                    intersectPoints_Z2 = [ intersectPoints_Z2 , Plane_Z2.PointA , Plane_Z2.PointB , Plane_Z2.PointC , Plane_Z2.PointD ];
                end
                intersectCoordTemp = [ intersectCoordTemp , intersectPoints_Z2 ];
            end
            
            intersectCoordTemp = [ intersectCoordTemp , Tetragon.PointA , Tetragon.PointB , Tetragon.PointC , Tetragon.PointD ];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Processing the intersection Points
            for nr = 1 : size( intersectCoordTemp , 2 )
                
                % Removing the Points that are not inside the matrix cell
                if ( ( XimCube(1) - intersectCoordTemp(1,nr) ) > Epsilon ) || ( ( intersectCoordTemp(1,nr) - XimCube(2) ) > Epsilon )
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                if ( ( YimCube(1) - intersectCoordTemp(2,nr) ) > Epsilon ) || ( ( intersectCoordTemp(2,nr) - YimCube(2) ) > Epsilon )
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                if ( ( ZimCube(1) - intersectCoordTemp(3,nr) ) > Epsilon ) || ( ( intersectCoordTemp(3,nr) - ZimCube(2) ) > Epsilon )
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                
                % Removing the Points that are not inside the fracture cell
                isInside = Tetragon.Is_Point_Inside_Tetragon(intersectCoordTemp(:,nr), Epsilon);
                if isInside == 0
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                
                % Removing the Points that are repeated
                subtracted = ( intersectCoordTemp - intersectCoordTemp(:,nr) );
                subtracted = sqrt( subtracted(1,:).^2 + subtracted(2,:).^2 + subtracted(3,:).^2 );
                if ~isempty( find( subtracted([1:nr-1 , nr+1:end]) <= Epsilon , 1 ) )
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                
                % Writing the intersection Points into a new variable
                if ~isnan(intersectCoordTemp(1,nr))
                    intersectCoordFinal = [ intersectCoordFinal , intersectCoordTemp(:,nr) ];
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Ordering the intersection Points (either clockwise or counter clockwise)
            intersectCoordTemp = intersectCoordFinal;
            if size(intersectCoordTemp,2) > 2
                distanceToZero = sqrt( sum( intersectCoordTemp.^2 , 1 ) );
                
                % Setting the 1st Point
                [ ~ , index ] = min( distanceToZero );
                intersectCoordFinal      = zeros ( size(intersectCoordTemp,1) , size(intersectCoordTemp,2) );
                intersectCoordFinal(:,1) = intersectCoordTemp(:,index);
                intersectCoordTemp(:,index) = [];
                
                % Setting the 2nd Point
                distanceToRef = sqrt( sum ( ( intersectCoordTemp - intersectCoordFinal(:,1) ).^2 , 1 ) );
                [ ~ , index ] = min( distanceToRef );
                intersectCoordFinal(:,2) = intersectCoordTemp(:,index);
                intersectCoordTemp(:,index) = [];
                
                % Setting the rest of the Points based on the maximum angle between two close segments
                for m = 2 : length( distanceToZero ) - 2
                    arcCos_Theta = ones( size( intersectCoordTemp , 2 ) , 1 );
                    for n = 1 : size( intersectCoordTemp , 2 )
                        arcCos_Theta(n) = dot( intersectCoordFinal(:,m) - intersectCoordFinal(:,m-1) , ...
                            intersectCoordFinal(:,m) - intersectCoordTemp (:,n  ) ) ...
                            / norm( intersectCoordFinal(:,m) - intersectCoordFinal(:,m-1) ) ...
                            / norm( intersectCoordFinal(:,m) - intersectCoordTemp (:,n  ) );
                    end
                    [ ~ , index ] = min( arcCos_Theta );
                    intersectCoordFinal(:,m+1)  = intersectCoordTemp(:,index);
                    intersectCoordTemp(:,index) = [];
                end
                intersectCoordFinal(:,end) = intersectCoordTemp;
            end
        end
    end
end