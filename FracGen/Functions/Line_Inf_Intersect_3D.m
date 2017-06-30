%%%%%%%%%%%%%%%%%%%%%%%%  Line_Inf_Intersect_3D  %%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2017-01-27
% Modified on: 2017-02-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lineGeostatus, lineIntersectPoint] = Line_Inf_Intersect_3D( point1 , vec1 , point2 , vec2 , almostZero )

lineGeostatus.haveIntersect = NaN;
lineGeostatus.areParallel   = NaN;
lineGeostatus.areCollinear  = NaN;
lineIntersectPoint = [];

%% Three main types of geostatus: 1- Parallel (either collinear or not), 2- intersected, 3- Skew
if norm( cross(vec1,vec2) ) < almostZero
    % The lines are either parallel or even collinear.
    lineGeostatus.areParallel   = 1;
    lineGeostatus.haveIntersect = 0;
    
    % checking if there is any distance between the prarallel lines
    parallelDistance = norm( cross( (point2 - point1) , vec1 ) ) / norm(vec1);
    if parallelDistance > almostZero
        % The lines are parallel, no intersection occurs.
        lineGeostatus.areCollinear  = 0;
        lineGeostatus.haveIntersect = 0;
    else
        % The lines are collinear, no intersection is considered.
        lineGeostatus.areCollinear  = 1;
        lineGeostatus.haveIntersect = 0;
    end
    
else
    % The lines are either intersected or skew.
    lineGeostatus.areParallel  = 0;
    lineGeostatus.areCollinear = 0; 
    t = det( [ (point2 - point1), vec2 , cross(vec1,vec2) ] ) / norm(cross(vec1,vec2))^2;
    s = det( [ (point2 - point1), vec1 , cross(vec1,vec2) ] ) / norm(cross(vec1,vec2))^2;
    
    if norm( (point1 + t*vec1) - (point2 + s*vec2) ) < almostZero
        % Intersection occurs
        lineGeostatus.haveIntersect  = 1;
        lineIntersectPoint = point1 + vec1.*t;
    else
        % The lines are skew
        lineGeostatus.haveIntersect  = 0;
    end
end

end