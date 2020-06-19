% Class of hexahedron for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim Reservoir Simulator
% Author: Mousa HosseiniMehr
% TU Delft
% Created: 2020-03-09
% Last modified: 2020-06-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef hexahedron_DARSim < polyhedron_DARSim
    properties
        NW_Top
        SW_Top
        SE_Top
        NE_Top
        NW_Bot
        SW_Bot
        SE_Bot
        NE_Bot
        Corners
        Centroid
        Face
        Volume
    end
    methods
        %%
        function obj = hexahedron_DARSim(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot)
            obj.NumOfVertex = 8;
            obj.NumOfFace = 6;
            obj.Face = tetragon_DARSim.empty;
            if nargin >= 8
                obj.InitializeHexahedron(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot);
            end
        end
        %%
        function InitializeHexahedron(obj,NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot)
            obj.AddCorners(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot);
            obj.InitializeFaces();
        end
        %%
        function AddCorners(obj,NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot)
            obj.NW_Top = NW_Top;
            obj.SW_Top = SW_Top;
            obj.SE_Top = SE_Top;
            obj.NE_Top = NE_Top;
            obj.NW_Bot = NW_Bot;
            obj.SW_Bot = SW_Bot;
            obj.SE_Bot = SE_Bot;
            obj.NE_Bot = NE_Bot;
            obj.Corners = [NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot];
        end
        %%
        function InitializeFaces(obj)
            obj.Face(1,1) = tetragon_DARSim(obj.NW_Bot,obj.NE_Bot,obj.NE_Top,obj.NW_Top,'North');
            obj.Face(2,1) = tetragon_DARSim(obj.SW_Bot,obj.SE_Bot,obj.SE_Top,obj.SW_Top,'South');
            obj.Face(3,1) = tetragon_DARSim(obj.SE_Bot,obj.NE_Bot,obj.NE_Top,obj.SE_Top,'East');
            obj.Face(4,1) = tetragon_DARSim(obj.SW_Bot,obj.NW_Bot,obj.NW_Top,obj.SW_Top,'West');
            obj.Face(5,1) = tetragon_DARSim(obj.SW_Top,obj.SE_Top,obj.NE_Top,obj.NW_Top,'Top');
            obj.Face(6,1) = tetragon_DARSim(obj.SW_Bot,obj.SE_Bot,obj.NE_Bot,obj.NW_Bot,'Bottom');
        end
        %%
        function AvgDistance = Obtain_Average_Distance_Plane_From_Hexahedron(obj, Plane, Refinement)
            % This functions calculates the average distance between a hexahedron and a plane.
            
            % PlaneEquation has 4 set of values of "a,b,c,d" and follows the equation
            % pattern "ax+by+cz=d".
            
            % "Refinement" is number of devisions in every direction. For example,
            % if the refinement number is 7, the cube is devided into 7^3 = 343 smaller
            % cubes to calculate the average distance.
            
            AvgDistance = 0;
            Refinement = round(max(Refinement, 2));
            
            a = Plane.Eq.a;
            b = Plane.Eq.b;
            c = Plane.Eq.c;
            d = Plane.Eq.d;
            
            Vec1_mini = ( obj.SE_Bot - obj.SW_Bot ) .* (0:Refinement) / Refinement;
            Vec2_mini = ( obj.NW_Bot - obj.SW_Bot ) .* (0:Refinement) / Refinement;
            Vec3_mini = ( obj.SW_Top - obj.SW_Bot ) .* (0:Refinement) / Refinement;

            for i = 1 : Refinement+1
                for j = 1 : Refinement+1
                    PointsInHexahedron = obj.SW_Bot + Vec1_mini(i) + Vec2_mini(j) + Vec3_mini;
                    AvgDistance = AvgDistance + sum ( abs( a * PointsInHexahedron(1,:) + b * PointsInHexahedron(2,:) + c * PointsInHexahedron(3,:) - d ) );
                end
            end

            AvgDistance = AvgDistance / sqrt( a^2 + b^2 + c^2 ) / (Refinement+1)^3 ;
        end
        %%
        function [intersectCoordFinal, intersectCoordTemp , areCoplanar] = Obtain_Hexahedron_Tetragon_Intersection(obj, Tetragon, Epsilon)
            % Initializing some variables
            doNotContinue       = 0;
            areCoplanar         = 0;
            intersectCoordTemp  = []; 
            intersectCoordFinal = [];
            
            for i= 1 : obj.NumOfFace
                if doNotContinue == 0
                    [Geostatus, intersectPoints] = Tetragon.Obtain_Tetragon_Tetragon_Intersection( obj.Face(i) , Epsilon );
                    if ( Geostatus.areCoplanar == 1 )
                        doNotContinue      = 1;
                        areCoplanar        = 1;
                        intersectPoints = [ intersectPoints , obj.Face(i).PointA , obj.Face(i).PointB , obj.Face(i).PointC , obj.Face(i).PointD ];
                    end
                    intersectCoordTemp = [ intersectCoordTemp , intersectPoints ];
                end
            end
            
            intersectCoordTemp = [ intersectCoordTemp , Tetragon.PointA , Tetragon.PointB , Tetragon.PointC , Tetragon.PointD ];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Processing the intersection Points
            for nr = 1 : size( intersectCoordTemp , 2 )
                
                % Removing the Points that are not inside the hexahedron
                if ( ( min(obj.Corners(1,:)) - intersectCoordTemp(1,nr) ) > Epsilon ) || ( ( intersectCoordTemp(1,nr) - max(obj.Corners(1,:)) ) > Epsilon )
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                if ( ( min(obj.Corners(2,:)) - intersectCoordTemp(2,nr) ) > Epsilon ) || ( ( intersectCoordTemp(2,nr) - max(obj.Corners(2,:)) ) > Epsilon )
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                if ( ( min(obj.Corners(3,:)) - intersectCoordTemp(3,nr) ) > Epsilon ) || ( ( intersectCoordTemp(3,nr) - max(obj.Corners(3,:)) ) > Epsilon )
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                
                % Removing the Points that are not inside the tetragon
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

