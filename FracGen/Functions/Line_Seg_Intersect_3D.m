%%%%%%%%%%%%%%%%%%%%%%%%  Line_Seg_Intersect_3D  %%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2017-01-02
% Modified on: 2017-02-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lineGeostatus, lineIntersectPoint] = Line_Seg_Intersect_3D(L1A,L1B,L2A,L2B,almostZero)

lineGeostatus.haveIntersect = NaN;
lineGeostatus.areParallel   = NaN;
lineGeostatus.areCollinear  = NaN;
lineIntersectPoint = [];

V1 = L1B - L1A;
V2 = L2B - L2A;

%% Three main types of geostatus: 1- Parallel (either collinear or not), 2- intersected(either with intersection points or not), 3- Skew
if ( norm( cross(V1,V2) ) / norm(V1) / norm(V2) ) < almostZero
    % The lines are either parallel or even collinear.
    lineGeostatus.areParallel   = 1;
    
    % checking if there is any distance between the prarallel line segment
    parallelDistance = norm( cross( (L2A - L1A) , V1 ) ) / norm(V1);
    if parallelDistance > almostZero
        % The lines are parallel, no intersection occurs.
        lineGeostatus.haveIntersect = 0;
        lineGeostatus.areCollinear  = 0;
        
    else
        % The lines are collinear.
        lineGeostatus.areCollinear = 1;
        intersectNr = 0;
        
        % Checking if point L2A is an intersection or not
        arcCos_Theta_L2A_1 = dot( L1B-L1A , L2A-L1A ) / ( norm(L1B-L1A) * norm(L2A-L1A) );
        arcCos_Theta_L2A_2 = dot( L1A-L1B , L2A-L1B ) / ( norm(L1A-L1B) * norm(L2A-L1B) );
        if ( arcCos_Theta_L2A_1 * arcCos_Theta_L2A_2 > almostZero ) && ( intersectNr <= 2 )
            intersectNr = intersectNr + 1;
            lineGeostatus.haveIntersect = 1;
            lineIntersectPoint = [ lineIntersectPoint , L2A ];
        end
        
        % Checking if point L2B is an intersection or not
        arcCos_Theta_L2B_1 = dot( L1B-L1A , L2B-L1A ) / ( norm(L1B-L1A) * norm(L2B-L1A) );
        arcCos_Theta_L2B_2 = dot( L1A-L1B , L2B-L1B ) / ( norm(L1A-L1B) * norm(L2B-L1B) );
        if ( arcCos_Theta_L2B_1 * arcCos_Theta_L2B_2 > almostZero ) && ( intersectNr <= 2 )
            intersectNr = intersectNr + 1;
            lineGeostatus.haveIntersect = 1;
            lineIntersectPoint = [ lineIntersectPoint , L2B ];
        end
        
        % Checking if point L1A is an intersection or not
        arcCos_Theta_L1A_1 = dot( L2B-L2A , L1A-L2A ) / ( norm(L2B-L2A) * norm(L1A-L2A) );
        arcCos_Theta_L1A_2 = dot( L2A-L2B , L1A-L2B ) / ( norm(L2A-L2B) * norm(L1A-L2B) );
        if ( arcCos_Theta_L1A_1 * arcCos_Theta_L1A_2 > almostZero ) && ( intersectNr <= 2 )
            intersectNr = intersectNr + 1;
            lineGeostatus.haveIntersect = 1;
            lineIntersectPoint = [ lineIntersectPoint , L1A ];
        end
        
        % Checking if point L1B is an intersection or not
        arcCos_Theta_L1B_1 = dot( L2B-L2A , L1B-L2A ) / ( norm(L2B-L2A) * norm(L1B-L2A) );
        arcCos_Theta_L1B_2 = dot( L2A-L2B , L1B-L2B ) / ( norm(L2A-L2B) * norm(L1B-L2B) );
        if ( arcCos_Theta_L1B_1 * arcCos_Theta_L1B_2 > almostZero ) && ( intersectNr <= 2 )
            intersectNr = intersectNr + 1;
            lineGeostatus.haveIntersect = 1;
            lineIntersectPoint = [ lineIntersectPoint , L1B ];
        end
        
    end
    
else
    % The lines are either intersected or skew.
    lineGeostatus.areParallel  = 0;
    lineGeostatus.areCollinear = 0;
    t = det( [ (L2A-L1A), V2 , cross(V1,V2) ] ) / norm(cross(V1,V2))^2;
    s = det( [ (L2A-L1A), V1 , cross(V1,V2) ] ) / norm(cross(V1,V2))^2;
    
    if norm( (L1A + t*V1) - (L2A + s*V2) ) < almostZero
        % Intersection occurs
        if ( t >= almostZero && t <= 1 + almostZero ) && ( s >= almostZero && s <= 1 + almostZero )
            % Intersection point lies in the line segments
            lineGeostatus.haveIntersect  = 1;
            lineIntersectPoint = L1A + V1.*t;
        else
            % Intersection point lies out of the line segments
            lineGeostatus.haveIntersect  = 0;
        end
    else
        % The lines are skew
        lineGeostatus.haveIntersect  = 0;
    end
end
            
end