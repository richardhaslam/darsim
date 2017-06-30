%%%%%%%%%%%%%%%%%%%  Is_Point_Inside_Rectangle_3D  %%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2017-01-02
% Modified on: 2016-02-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isInside = Is_Point_Inside_Rectangle_3D( Plane , Point , almostZero )
%%
isInside  = NaN;

Plane.PointM = ( Plane.PointB + Plane.PointD )/2;
almostZero = almostZero * 1e3;

if dot( cross((Plane.PointA-Plane.PointM),(Plane.PointB-Plane.PointM)) , (Point-Plane.PointM) ) > almostZero
    % The Point is not on the plane
    isInside = 0;  return;
end

[lineGeostatus, lineIntersectPoint] = Line_Seg_Intersect_3D( Plane.PointA , Plane.PointB , Plane.PointM , Point , almostZero );
if     ( lineGeostatus.haveIntersect == 0 )                ,  isInside = 1;
elseif ( norm (lineIntersectPoint - Point ) < almostZero ) ,  isInside = 1;
else                                                       ,  isInside = 0;  return;  end
    
[lineGeostatus, lineIntersectPoint] = Line_Seg_Intersect_3D( Plane.PointB , Plane.PointC , Plane.PointM , Point , almostZero );
if     ( lineGeostatus.haveIntersect == 0 )                ,  isInside = 1;
elseif ( norm (lineIntersectPoint - Point ) < almostZero ) ,  isInside = 1;
else                                                       ,  isInside = 0;  return;  end

[lineGeostatus, lineIntersectPoint] = Line_Seg_Intersect_3D( Plane.PointC , Plane.PointD , Plane.PointM , Point , almostZero );
if     ( lineGeostatus.haveIntersect == 0 )                ,  isInside = 1;
elseif ( norm (lineIntersectPoint - Point ) < almostZero ) ,  isInside = 1;
else                                                       ,  isInside = 0;  return;  end

[lineGeostatus, lineIntersectPoint] = Line_Seg_Intersect_3D( Plane.PointD , Plane.PointA , Plane.PointM , Point , almostZero );
if     ( lineGeostatus.haveIntersect == 0 )                ,  isInside = 1;
elseif ( norm (lineIntersectPoint - Point ) < almostZero ) ,  isInside = 1;
else                                                       ,  isInside = 0;  return;  end
            
end