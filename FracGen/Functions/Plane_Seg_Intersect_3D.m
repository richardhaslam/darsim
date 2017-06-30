%%%%%%%%%%%%%%%%%%%%%%%  Plane_Seg_Intersect_3D  %%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2017-01-02
% Modified on: 2017-02-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [planeGeostatus, planeIntersectPoints] = Plane_Seg_Intersect_3D(Plane1,Plane2,almostZero)

planeGeostatus.areParallel   = NaN;
planeGeostatus.areCoplanar   = NaN;
planeGeostatus.haveIntersect = NaN;
planeIntersectPoints         = [];

% AB, AD and AC of both plane segments
Plane1.AB_vec = Plane1.PointB - Plane1.PointA;  Plane1.AD_vec = Plane1.PointD - Plane1.PointA;  Plane1.AC_vec = Plane1.PointC - Plane1.PointA;
Plane2.AB_vec = Plane2.PointB - Plane2.PointA;  Plane2.AD_vec = Plane2.PointD - Plane2.PointA;  Plane2.AC_vec = Plane2.PointC - Plane2.PointA;

% Normal vector of both plane segments
Plane1.n_vec = cross( Plane1.AB_vec , Plane1.AD_vec );  Plane1.n_vec = Plane1.n_vec / norm( Plane1.n_vec );
Plane2.n_vec = cross( Plane2.AB_vec , Plane2.AD_vec );  Plane2.n_vec = Plane2.n_vec / norm( Plane2.n_vec );

% RHS of plane equation for both plane segments correspondent with ax+by+cz=d
Plane1.Equation.d = dot( Plane1.PointA , Plane1.n_vec );
Plane2.Equation.d = dot( Plane2.PointA , Plane2.n_vec );

%% Two main types of geostatus: 1- Parallel (either coplanar or not), 2- intersected(either with intersection Points or not)
if norm( cross( Plane1.n_vec , Plane2.n_vec ) ) < almostZero
    % The planes are parallel or even coplanar:
    
    % checking if Point Plane2.PointA lies on the 2nd plane (coplanar) or not (only parallel)
    if abs( dot( Plane2.PointA , Plane1.n_vec ) - Plane1.Equation.d ) > almostZero
        % The planes are parallel. No intersection is considered.
        planeGeostatus.areParallel   = 1;
        planeGeostatus.areCoplanar   = 0;
        planeGeostatus.haveIntersect = 0;

    else
        % The planes are coplanar.
        planeGeostatus.areParallel   = 1;
        planeGeostatus.areCoplanar   = 1;
        
        Line_1st = [ Plane1.PointA , Plane1.PointB , Plane1.PointC , Plane1.PointD
                     Plane1.PointB , Plane1.PointC , Plane1.PointD , Plane1.PointA ];
        Line_2nd = [ Plane2.PointA , Plane2.PointB , Plane2.PointC , Plane2.PointD
                     Plane2.PointB , Plane2.PointC , Plane2.PointD , Plane2.PointA ];
        
        % Checking the intersections between each two edges of plane segments         
        for i = 1 : size(Line_1st,2)
            intersectNr = 0;
            for j = 1 : size(Line_2nd,2)
                [lineGeostatus, lineIntersectPoint] = Line_Seg_Intersect_3D( Line_1st(1:3,i) , Line_1st(4:6,i) , Line_2nd(1:3,j) , Line_2nd(4:6,j) , almostZero );
                if lineGeostatus.haveIntersect == 1
                    planeGeostatus.haveIntersect = 1;
                    intersectNr = intersectNr + 1;
                    planeIntersectPoints = [planeIntersectPoints , lineIntersectPoint];
                end
                if intersectNr == 2,  continue;  end  % each edge of 1st plane segment can have intersection with maximum two edges of the 2nd plane segment
            end
        end
           
    end
  
else % The planes will intersect.

    planeGeostatus.areParallel = 0;
    planeGeostatus.areCoplanar = 0;
        
    % Obtaining the unit vector of the intersection line
    intL_Vec = cross( Plane1.n_vec , Plane2.n_vec );
    intL_Vec = intL_Vec / norm(intL_Vec);
    
    % Obtaining a Point on the intersection line with Z=0 if possible
    P2M = (Plane2.PointA + Plane2.PointC )/2;
    Axy = [ Plane1.n_vec(1) , Plane1.n_vec(2) ; Plane2.n_vec(1)   , Plane2.n_vec(2) ];
    RHS = [ Plane1.Equation.d - Plane1.n_vec(3)*P2M(3)
            Plane2.Equation.d - Plane2.n_vec(3)*P2M(3) ];
    if abs(det(Axy)) > almostZero
        intL_Point0    = zeros(3,1);
        Unknowns       = Axy \ RHS;
        intL_Point0(1) = Unknowns(1);
        intL_Point0(2) = Unknowns(2);
        intL_Point0(3) = P2M(3);
    else
        % Or, obtaining a Point on the intersection line with Y=0 if possible
        Axz = [ Plane1.n_vec(1) , Plane1.n_vec(3) ; Plane2.n_vec(1)  , Plane2.n_vec(3) ];
        RHS = [ Plane1.Equation.d - Plane1.n_vec(2)*P2M(2)
                Plane2.Equation.d - Plane2.n_vec(2)*P2M(2) ];
        if abs(det(Axz)) > almostZero
            intL_Point0    = zeros(3,1);
            Unknowns       = Axz \ RHS;
            intL_Point0(1) = Unknowns(1);
            intL_Point0(3) = Unknowns(2);
            intL_Point0(2) = P2M(2);
        else
            % Or, obtaining a Point on the intersection line with X=0 if possible
            Ayz = [ Plane1.n_vec(2) , Plane1.n_vec(3) ; Plane2.n_vec(2)  , Plane2.n_vec(3) ];
            RHS = [ Plane1.Equation.d - Plane1.n_vec(1)*P2M(1)
                    Plane2.Equation.d - Plane2.n_vec(1)*P2M(1) ];
            intL_Point0    = zeros(3,1);
            Unknowns       = Ayz \ RHS;
            intL_Point0(2) = Unknowns(1);
            intL_Point0(3) = Unknowns(2);
            intL_Point0(1) = P2M(1);
        end
    end

    % Assuming two end-Points for the intersection line
    intL_A = intL_Point0 - intL_Vec * max( norm(Plane1.AC_vec) , norm(Plane2.AC_vec) ) * 100;
    intL_B = intL_Point0 + intL_Vec * max( norm(Plane1.AC_vec) , norm(Plane2.AC_vec) ) * 100;
    
    % Checking the intersection between the obtained line segment and each side of each plane segment
    Point_1st = [ Plane1.PointA , Plane1.PointB , Plane1.PointC , Plane1.PointD , Plane2.PointA , Plane2.PointB , Plane2.PointC , Plane2.PointD ];
    Point_2nd = [ Plane1.PointB , Plane1.PointC , Plane1.PointD , Plane1.PointA , Plane2.PointB , Plane2.PointC , Plane2.PointD , Plane2.PointA ];
    
    intersectNr = 0;
    for i = 1 : size(Point_1st,2)
        [lineGeostatus, lineIntersectPoint] = Line_Seg_Intersect_3D( Point_1st(:,i) , Point_2nd(:,i) , intL_A , intL_B , almostZero );
        if lineGeostatus.haveIntersect == 1
            planeGeostatus.haveIntersect = 1;
            intersectNr = intersectNr + 1;
            planeIntersectPoints = [planeIntersectPoints , lineIntersectPoint];
        end
        if intersectNr == 4,  continue;  end % intersection line can have intersection with maximum two edges of each plane segment (four intersections max)
    end
   
end

end