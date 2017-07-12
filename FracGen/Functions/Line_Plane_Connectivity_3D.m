%%%%%%%%%%%%%%%%%%%%%%  Frac_Frac_Connectivity_3D  %%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2017-01-27
% Modified on: 2016-01-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Af , aveDist , Collinearity] = Line_Plane_Connectivity_3D( Plane , Point1 , Point2 , almostZero )

Plane.AB = Plane.PointB - Plane.PointA;
Plane.BC = Plane.PointC - Plane.PointB;
Plane.CD = Plane.PointD - Plane.PointC;
Plane.DA = Plane.PointA - Plane.PointD;

vec12 = Point2 - Point1;
Collinearity = 0;
%% Calculating Af of the line segment inside the plane segment
Af = 2 * norm ( Point1 - Point2 );

[lineGeostatus_AB, lineIntersectPoint_AB] = Line_Inf_Intersect_3D( Plane.PointA , Plane.AB , Point1 , vec12 , almostZero );
[lineGeostatus_BC, lineIntersectPoint_BC] = Line_Inf_Intersect_3D( Plane.PointB , Plane.BC , Point1 , vec12 , almostZero );
[lineGeostatus_CD, lineIntersectPoint_CD] = Line_Inf_Intersect_3D( Plane.PointC , Plane.CD , Point1 , vec12 , almostZero );
[lineGeostatus_DA, lineIntersectPoint_DA] = Line_Inf_Intersect_3D( Plane.PointD , Plane.DA , Point1 , vec12 , almostZero );

if ( lineGeostatus_AB.areCollinear == 1 ) || ( lineGeostatus_BC.areCollinear == 1 ) || ...
   ( lineGeostatus_CD.areCollinear == 1 ) || ( lineGeostatus_DA.areCollinear == 1 )
    Collinearity = 1;
    Af = Af /2;
end

%% Obtaining aveDist of the plane segment from the line segment

% Parallel situation
if ( lineGeostatus_AB.areParallel == 1 ) || ( lineGeostatus_CD.areParallel == 1 )
    % The line segment is parallel to AB and CD
    Dist1 = norm( cross( (Point1 - Plane.PointA) , vec12 ) ) / norm(vec12);
    Dist2 = norm( cross( (Point1 - Plane.PointC) , vec12 ) ) / norm(vec12);
    if abs(Dist1+Dist2-norm(Plane.DA)) > almostZero,  error('The summation of distances of the parallel line from AB and CD is not correct!');  end
    aveDist = ( Dist1^2 + Dist2^2 ) / ( 2*Dist1 + 2*Dist2 );
    return;
end
if ( lineGeostatus_BC.areParallel == 1 ) || ( lineGeostatus_DA.areParallel == 1 )
    % The line segment is parallel to BC and DA
    Dist1 = norm( cross( (Point1 - Plane.PointB) , vec12 ) ) / norm(vec12);
    Dist2 = norm( cross( (Point1 - Plane.PointD) , vec12 ) ) / norm(vec12);
    if abs(Dist1+Dist2-norm(Plane.AB)) > almostZero,  error('The summation of distances of the parallel line from BC and DA is not correct!');  end
    aveDist = ( Dist1^2 + Dist2^2 ) / ( 2*Dist1 + 2*Dist2 );
    return;
end

% Non-parallel situation
isInside_AB = Is_Point_Inside_Rectangle_3D( Plane , lineIntersectPoint_AB , almostZero );
isInside_BC = Is_Point_Inside_Rectangle_3D( Plane , lineIntersectPoint_BC , almostZero );
isInside_CD = Is_Point_Inside_Rectangle_3D( Plane , lineIntersectPoint_CD , almostZero );
isInside_DA = Is_Point_Inside_Rectangle_3D( Plane , lineIntersectPoint_DA , almostZero );

if sum( [ isInside_AB , isInside_BC , isInside_CD , isInside_DA ] ) == 4
    % The line segment is on diagonal of plane segment
    aveDist = norm(Plane.AB) * norm(Plane.BC) / ( 3 * sqrt( norm(Plane.AB)^2 + norm(Plane.BC)^2 ) );
    return;
end

if ( sum( [ isInside_AB , isInside_CD ] ) == 2 )
    % The extension of the line segment intersects with two parallel sides of the plane segment (AB and CD)
    if norm( Plane.PointA - lineIntersectPoint_DA ) < norm( Plane.PointD - lineIntersectPoint_DA )
        Area1 = Triangle_Area_3D( lineIntersectPoint_DA , Plane.PointD , lineIntersectPoint_CD );  % Large Triangle
        Area2 = Triangle_Area_3D( lineIntersectPoint_DA , Plane.PointA , lineIntersectPoint_AB );  % Small Triangle
        Area3 = Triangle_Area_3D( lineIntersectPoint_BC , Plane.PointB , lineIntersectPoint_AB );  % Large Triangle
        Area4 = Triangle_Area_3D( lineIntersectPoint_BC , Plane.PointC , lineIntersectPoint_CD );  % Small Triangle
        Length1 = norm( lineIntersectPoint_CD - Plane.PointD );   Width1 = norm( Plane.PointD - lineIntersectPoint_CD );  % Large Triangle
        Length2 = norm( lineIntersectPoint_AB - Plane.PointA );   Width2 = norm( Plane.PointA - lineIntersectPoint_AB );  % Small Triangle
        Length3 = norm( lineIntersectPoint_AB - Plane.PointB );   Width3 = norm( Plane.PointB - lineIntersectPoint_AB );  % Large Triangle
        Length4 = norm( lineIntersectPoint_CD - Plane.PointC );   Width4 = norm( Plane.PointC - lineIntersectPoint_CD );  % Small Triangle
    else
        Area1 = Triangle_Area_3D( lineIntersectPoint_DA , Plane.PointA , lineIntersectPoint_AB );  % Large Triangle
        Area2 = Triangle_Area_3D( lineIntersectPoint_DA , Plane.PointD , lineIntersectPoint_CD );  % Small Triangle
        Area3 = Triangle_Area_3D( lineIntersectPoint_BC , Plane.PointC , lineIntersectPoint_CD );  % Large Triangle
        Area4 = Triangle_Area_3D( lineIntersectPoint_BC , Plane.PointB , lineIntersectPoint_AB );  % Small Triangle
        Length1 = norm( lineIntersectPoint_AB - Plane.PointA );   Width1 = norm( Plane.PointA - lineIntersectPoint_AB );  % Large Triangle
        Length2 = norm( lineIntersectPoint_CD - Plane.PointD );   Width2 = norm( Plane.PointD - lineIntersectPoint_CD );  % Small Triangle
        Length3 = norm( lineIntersectPoint_CD - Plane.PointC );   Width3 = norm( Plane.PointC - lineIntersectPoint_CD );  % Large Triangle
        Length4 = norm( lineIntersectPoint_AB - Plane.PointB );   Width4 = norm( Plane.PointB - lineIntersectPoint_AB );  % Small Triangle
    end
    Dist1 = Length1 * Width1 / ( 3 * sqrt( Length1^2 + Width1^2 ) );  if Area1==0,  Dist1=0;  end
    Dist2 = Length2 * Width2 / ( 3 * sqrt( Length2^2 + Width2^2 ) );  if Area2==0,  Dist2=0;  end
    Dist3 = Length3 * Width3 / ( 3 * sqrt( Length3^2 + Width3^2 ) );  if Area3==0,  Dist3=0;  end
    Dist4 = Length4 * Width4 / ( 3 * sqrt( Length4^2 + Width4^2 ) );  if Area4==0,  Dist4=0;  end
    aveDist = ( Area1*Dist1 - Area2*Dist2 + Area3*Dist3 - Area4*Dist4 ) / ( Area1 - Area2 + Area3 - Area4 );
    return;
end

if ( sum( [ isInside_BC , isInside_DA ] ) == 2 )
    % The extension of the line segment intersects with two parallel sides of the plane segment (BC and DA)
    if norm( Plane.PointA - lineIntersectPoint_AB ) < norm( Plane.PointB - lineIntersectPoint_AB )
        Area1 = Triangle_Area_3D( lineIntersectPoint_AB , Plane.PointB , lineIntersectPoint_BC );  % Large Triangle
        Area2 = Triangle_Area_3D( lineIntersectPoint_AB , Plane.PointA , lineIntersectPoint_DA );  % Small Triangle
        Area3 = Triangle_Area_3D( lineIntersectPoint_CD , Plane.PointD , lineIntersectPoint_DA );  % Large Triangle
        Area4 = Triangle_Area_3D( lineIntersectPoint_CD , Plane.PointC , lineIntersectPoint_BC );  % Small Triangle
        Length1 = norm( lineIntersectPoint_AB - Plane.PointB );   Width1 = norm( Plane.PointB - lineIntersectPoint_BC );  % Large Triangle
        Length2 = norm( lineIntersectPoint_AB - Plane.PointA );   Width2 = norm( Plane.PointA - lineIntersectPoint_DA );  % Small Triangle
        Length3 = norm( lineIntersectPoint_CD - Plane.PointD );   Width3 = norm( Plane.PointD - lineIntersectPoint_DA );  % Large Triangle
        Length4 = norm( lineIntersectPoint_CD - Plane.PointC );   Width4 = norm( Plane.PointC - lineIntersectPoint_BC );  % Small Triangle
    else
        Area1 = Triangle_Area_3D( lineIntersectPoint_AB , Plane.PointA , lineIntersectPoint_DA );  % Large Triangle
        Area2 = Triangle_Area_3D( lineIntersectPoint_AB , Plane.PointB , lineIntersectPoint_BC );  % Small Triangle
        Area3 = Triangle_Area_3D( lineIntersectPoint_CD , Plane.PointC , lineIntersectPoint_BC );  % Large Triangle
        Area4 = Triangle_Area_3D( lineIntersectPoint_CD , Plane.PointD , lineIntersectPoint_DA );  % Small Triangle
        Length1 = norm( lineIntersectPoint_AB - Plane.PointA );   Width1 = norm( Plane.PointA - lineIntersectPoint_DA );  % Large Triangle
        Length2 = norm( lineIntersectPoint_AB - Plane.PointB );   Width2 = norm( Plane.PointB - lineIntersectPoint_BC );  % Small Triangle
        Length3 = norm( lineIntersectPoint_CD - Plane.PointC );   Width3 = norm( Plane.PointC - lineIntersectPoint_BC );  % Large Triangle
        Length4 = norm( lineIntersectPoint_CD - Plane.PointD );   Width4 = norm( Plane.PointD - lineIntersectPoint_DA );  % Small Triangle
    end
    Dist1 = Length1 * Width1 / ( 3 * sqrt( Length1^2 + Width1^2 ) );  if Area1==0,  Dist1=0;  end
    Dist2 = Length2 * Width2 / ( 3 * sqrt( Length2^2 + Width2^2 ) );  if Area2==0,  Dist2=0;  end
    Dist3 = Length3 * Width3 / ( 3 * sqrt( Length3^2 + Width3^2 ) );  if Area3==0,  Dist3=0;  end
    Dist4 = Length4 * Width4 / ( 3 * sqrt( Length4^2 + Width4^2 ) );  if Area4==0,  Dist4=0;  end
    aveDist = ( Area1*Dist1 - Area2*Dist2 + Area3*Dist3 - Area4*Dist4 ) / ( Area1 - Area2 + Area3 - Area4 );
    return;
end

if ( sum( [ isInside_AB , isInside_BC ] ) == 2 )
    % The extension of the line segment intersects with two cornering sides of the plane segment (AB and BC)
    Area1 = Triangle_Area_3D( lineIntersectPoint_CD , Plane.PointD , lineIntersectPoint_DA );  % Largest Triangle
    Area2 = Triangle_Area_3D( lineIntersectPoint_AB , Plane.PointB , lineIntersectPoint_BC );  % Small Triangle Inside
    Area3 = Triangle_Area_3D( lineIntersectPoint_DA , Plane.PointA , lineIntersectPoint_AB );  % Outer Triangle 1
    Area4 = Triangle_Area_3D( lineIntersectPoint_BC , Plane.PointC , lineIntersectPoint_CD );  % Outer Triangle 2
    Length1 = norm( lineIntersectPoint_CD - Plane.PointD );   Width1 = norm( Plane.PointD - lineIntersectPoint_DA );  % Large Triangle
    Length2 = norm( lineIntersectPoint_AB - Plane.PointB );   Width2 = norm( Plane.PointB - lineIntersectPoint_BC );  % Small Triangle
    Length3 = norm( lineIntersectPoint_DA - Plane.PointA );   Width3 = norm( Plane.PointA - lineIntersectPoint_AB );  % Large Triangle
    Length4 = norm( lineIntersectPoint_BC - Plane.PointC );   Width4 = norm( Plane.PointC - lineIntersectPoint_CD );  % Small Triangle
    
elseif ( sum( [ isInside_BC , isInside_CD ] ) == 2 )
    % The extension of the line segment intersects with two cornering sides of the plane segment (BC and CD)
    Area1  = Triangle_Area_3D( lineIntersectPoint_DA , Plane.PointA , lineIntersectPoint_AB );  % Largest Triangle
    Area2 = Triangle_Area_3D( lineIntersectPoint_BC , Plane.PointC , lineIntersectPoint_CD );  % Small Triangle Inside
    Area3 = Triangle_Area_3D( lineIntersectPoint_AB , Plane.PointB , lineIntersectPoint_BC );  % Outer Triangle 1
    Area4 = Triangle_Area_3D( lineIntersectPoint_CD , Plane.PointD , lineIntersectPoint_DA );  % Outer Triangle 2
    Length1 = norm( lineIntersectPoint_DA - Plane.PointA );   Width1 = norm( Plane.PointA - lineIntersectPoint_AB );  % Large Triangle
    Length2 = norm( lineIntersectPoint_BC - Plane.PointC );   Width2 = norm( Plane.PointC - lineIntersectPoint_CD );  % Small Triangle
    Length3 = norm( lineIntersectPoint_AB - Plane.PointB );   Width3 = norm( Plane.PointB - lineIntersectPoint_BC );  % Large Triangle
    Length4 = norm( lineIntersectPoint_CD - Plane.PointD );   Width4 = norm( Plane.PointD - lineIntersectPoint_DA );  % Small Triangle
    
elseif ( sum( [ isInside_CD , isInside_DA ] ) == 2 )
    % The extension of the line segment intersects with two cornering sides of the plane segment (CD and DA)
    Area1 = Triangle_Area_3D( lineIntersectPoint_AB , Plane.PointB , lineIntersectPoint_BC );  % Largest Triangle
    Area2 = Triangle_Area_3D( lineIntersectPoint_CD , Plane.PointD , lineIntersectPoint_DA );  % Small Triangle Inside
    Area3 = Triangle_Area_3D( lineIntersectPoint_DA , Plane.PointA , lineIntersectPoint_AB );  % Outer Triangle 1
    Area4 = Triangle_Area_3D( lineIntersectPoint_BC , Plane.PointC , lineIntersectPoint_CD );  % Outer Triangle 2
    Length1 = norm( lineIntersectPoint_AB - Plane.PointB );   Width1 = norm( Plane.PointB - lineIntersectPoint_BC );  % Large Triangle
    Length2 = norm( lineIntersectPoint_CD - Plane.PointD );   Width2 = norm( Plane.PointD - lineIntersectPoint_DA );  % Small Triangle
    Length3 = norm( lineIntersectPoint_DA - Plane.PointA );   Width3 = norm( Plane.PointA - lineIntersectPoint_AB );  % Large Triangle
    Length4 = norm( lineIntersectPoint_BC - Plane.PointC );   Width4 = norm( Plane.PointC - lineIntersectPoint_CD );  % Small Triangle
    
elseif ( sum( [ isInside_DA , isInside_AB ] ) == 2 )
    % The extension of the line segment intersects with two cornering sides of the plane segment (DA and AB)
    Area1 = Triangle_Area_3D( lineIntersectPoint_BC , Plane.PointC , lineIntersectPoint_CD );  % Largest Triangle
    Area2 = Triangle_Area_3D( lineIntersectPoint_DA , Plane.PointA , lineIntersectPoint_AB );  % Small Triangle Inside
    Area3 = Triangle_Area_3D( lineIntersectPoint_AB , Plane.PointB , lineIntersectPoint_BC );  % Outer Triangle 1
    Area4 = Triangle_Area_3D( lineIntersectPoint_CD , Plane.PointD , lineIntersectPoint_DA );  % Outer Triangle 2
    Length1 = norm( lineIntersectPoint_BC - Plane.PointC );   Width1 = norm( Plane.PointC - lineIntersectPoint_CD );  % Large Triangle
    Length2 = norm( lineIntersectPoint_DA - Plane.PointA );   Width2 = norm( Plane.PointA - lineIntersectPoint_AB );  % Small Triangle
    Length3 = norm( lineIntersectPoint_AB - Plane.PointB );   Width3 = norm( Plane.PointB - lineIntersectPoint_BC );  % Large Triangle
    Length4 = norm( lineIntersectPoint_CD - Plane.PointD );   Width4 = norm( Plane.PointD - lineIntersectPoint_DA );  % Small Triangle
end

Dist1 = Length1 * Width1 / ( 3 * sqrt( Length1^2 + Width1^2 ) );  if Area1==0,  Dist1=0;  end
Dist2 = Length2 * Width2 / ( 3 * sqrt( Length2^2 + Width2^2 ) );  if Area2==0,  Dist2=0;  end
Dist3 = Length3 * Width3 / ( 3 * sqrt( Length3^2 + Width3^2 ) );  if Area3==0,  Dist3=0;  end
Dist4 = Length4 * Width4 / ( 3 * sqrt( Length4^2 + Width4^2 ) );  if Area4==0,  Dist4=0;  end
aveDist = ( Area1*Dist1 + Area2*Dist2 - Area3*Dist3 - Area4*Dist4 ) / ( Area1 + Area2 - Area3 - Area4 );

end