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
    end
end

