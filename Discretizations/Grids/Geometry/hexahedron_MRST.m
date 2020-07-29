% Class of hexahedron for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim Reservoir Simulator
% Author: Mousa HosseiniMehr
% TU Delft
% Created: 2020-03-09
% Last modified: 2020-06-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef hexahedron_MRST < hexahedron_DARSim
    properties
    end
    methods
        %%
        function obj = hexahedron_DARSim(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot)
            obj.NumOfVertex = 8;
            obj.Face = polygon_DARSim.empty;
            if nargin >= 8
                obj.AddVertices(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot);
            end
        end
        %%
        function AddFaces(obj,face)
            obj.NumOfFace = length(face);
            obj.Face = face;
        end
        %%
        function [intersectCoordFinal, intersectCoordTemp , areCoplanar] = Obtain_Hexahedron_Polygon_Intersection(obj, Polygon, Epsilon)
            % Initializing some variables
            doNotContinue       = 0;
            areCoplanar         = 0;
            intersectCoordTemp  = []; 
            intersectCoordFinal = [];
            
            for i= 1 : obj.NumOfFace
                if doNotContinue == 0
                    [Geostatus, intersectPoints] = Polygon.Obtain_Polygon_Polygon_Intersection( obj.Face(i) , Epsilon );
                    if ( Geostatus.areCoplanar == 1 )
                        doNotContinue      = 1;
                        areCoplanar        = 1;
                        intersectPoints = [ intersectPoints ; obj.Vertex ];
                    end
                    intersectCoordTemp = [ intersectCoordTemp ; intersectPoints ];
                end
            end
            
            intersectCoordTemp = [ intersectCoordTemp ; Polygon.Vertex ];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Processing the intersection Points
            for nr = 1 : size( intersectCoordTemp , 1 )
                
                % Removing the Points that are not inside the hexahedron
                if ( ( min(obj.Vertex(:,1)) - intersectCoordTemp(nr,1) ) > Epsilon ) || ( ( intersectCoordTemp(nr,1) - max(obj.Vertex(:,1)) ) > Epsilon )
                    intersectCoordTemp(nr,:) = [ NaN , NaN , NaN ];
                    continue;
                end
                if ( ( min(obj.Vertex(:,2)) - intersectCoordTemp(nr,2) ) > Epsilon ) || ( ( intersectCoordTemp(nr,2) - max(obj.Vertex(:,2)) ) > Epsilon )
                    intersectCoordTemp(nr,:) = [ NaN , NaN , NaN ];
                    continue;
                end
                if ( ( min(obj.Vertex(:,3)) - intersectCoordTemp(nr,3) ) > Epsilon ) || ( ( intersectCoordTemp(nr,3) - max(obj.Vertex(:,3)) ) > Epsilon )
                    intersectCoordTemp(nr,:) = [ NaN , NaN , NaN ];
                    continue;
                end
                
                % Removing the Points that are not inside the tetragon
                isInside = Polygon.Is_Point_Inside_Polygon(intersectCoordTemp(nr,:), Epsilon);
                if isInside == 0
                    intersectCoordTemp(nr,:) = [ NaN , NaN , NaN ];
                    continue;
                end
                
                % Removing the Points that are repeated
                subtracted = ( intersectCoordTemp - intersectCoordTemp(nr,:) );
                subtracted = sqrt( subtracted(:,1).^2 + subtracted(:,2).^2 + subtracted(:,3).^2 );
                if ~isempty( find( subtracted([1:nr-1 , nr+1:end]) <= Epsilon , 1 ) )
                    intersectCoordTemp(nr,:) = [ NaN , NaN , NaN ];
                    continue;
                end
                
                % Writing the intersection Points into a new variable
                if ~isnan(intersectCoordTemp(nr,1))
                    intersectCoordFinal = [ intersectCoordFinal ; intersectCoordTemp(nr,:) ];
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Ordering the intersection Points (either clockwise or counter clockwise)
            intersectCoordTemp = intersectCoordFinal;
            if size(intersectCoordTemp,1) > 2
                distanceToZero = sqrt( sum( intersectCoordTemp.^2 , 2 ) );
                
                % Setting the 1st Point
                [ ~ , index ] = min( distanceToZero );
                intersectCoordFinal      = zeros ( size(intersectCoordTemp) );
                intersectCoordFinal(1,:) = intersectCoordTemp(index,:);
                intersectCoordTemp(index,:) = [];
                
                % Setting the 2nd Point
                distanceToRef = sqrt( sum ( ( intersectCoordTemp - intersectCoordFinal(1,:) ).^2 , 2 ) );
                [ ~ , index ] = min( distanceToRef );
                intersectCoordFinal(2,:) = intersectCoordTemp(index,:);
                intersectCoordTemp(index,:) = [];
                
                % Setting the rest of the Points based on the maximum angle between two close segments
                for m = 2 : length( distanceToZero ) - 2
                    arcCos_Theta = ones( size( intersectCoordTemp , 1 ) , 1 );
                    for n = 1 : size( intersectCoordTemp , 1 )
                        arcCos_Theta(n) = dot( intersectCoordFinal(m,:) - intersectCoordFinal(m-1,:) , ...
                                               intersectCoordFinal(m,:) - intersectCoordTemp (n  ,:) ) ...
                            / norm( intersectCoordFinal(m,:) - intersectCoordFinal(m-1,:) ) ...
                            / norm( intersectCoordFinal(m,:) - intersectCoordTemp (n  ,:) );
                    end
                    [ ~ , index ] = min( arcCos_Theta );
                    intersectCoordFinal(m+1,:)  = intersectCoordTemp(index,:);
                    intersectCoordTemp(index,:) = [];
                end
                intersectCoordFinal(end,:) = intersectCoordTemp;
            end
        end
    end
end

