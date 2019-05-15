%%%%%%%%%%%%%%%%%%%  Average_Distance_Cube_From_Plane  %%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2017-01-07
% Modified on: 2017-01-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AvgDist = Average_Distance_Cube_From_Plane( PlaneEquation, CubeCornerStart , CubeCornerEnd , Refinement )

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

AvgDist = 0;

XimCube = linspace( CubeCornerStart(1), CubeCornerEnd(1), Refinement+1 );
YimCube = linspace( CubeCornerStart(2), CubeCornerEnd(2), Refinement+1 );
ZimCube = linspace( CubeCornerStart(3), CubeCornerEnd(3), Refinement+1 );

a= PlaneEquation.a; b= PlaneEquation.b; c= PlaneEquation.c; d= PlaneEquation.d;

for k = 1 : Refinement+1
    for j = 1 : Refinement+1
        AvgDist = AvgDist + sum ( abs( a * XimCube + b * YimCube(j) + c * ZimCube(k) - d ) );
    end
end

AvgDist = AvgDist / sqrt( a^2 + b^2 + c^2 ) / (Refinement+1)^3 ;

end