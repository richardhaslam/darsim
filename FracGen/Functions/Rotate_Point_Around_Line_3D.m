%%%%%%%%%%%%%%%%%%%%%  Rotate_Point_Around_Line_3D  %%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2017-01-07
% Modified on: 2017-01-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using the method in http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PointAfterRot = Rotate_Point_Around_Line_3D( PointPreRot, rotAxis , pointInAxis , rotAngle )

rotAxis = rotAxis / norm( rotAxis );
Theta   = rotAngle;

a = pointInAxis(1); u = rotAxis(1);
b = pointInAxis(2); v = rotAxis(2);
c = pointInAxis(3); w = rotAxis(3);
            
rotMat = [  u^2 + (v^2+w^2)*cos(Theta)         ,  u*v*(1-cos(Theta)) - w*sin(Theta)  ,  u*w*(1-cos(Theta)) + v*sin(Theta)  ,  ( a*(v^2+w^2) - u*(b*v+c*w) ) * (1-cos(Theta)) + (b*w-c*v)*sin(Theta)
            u*v*(1-cos(Theta)) + w*sin(Theta)  ,  v^2 + (u^2+w^2)*cos(Theta)         ,  v*w*(1-cos(Theta)) - u*sin(Theta)  ,  ( b*(u^2+w^2) - v*(a*u+c*w) ) * (1-cos(Theta)) + (c*u-a*w)*sin(Theta)
            u*w*(1-cos(Theta)) - v*sin(Theta)  ,  v*w*(1-cos(Theta)) + u*sin(Theta)  ,  w^2 + (u^2+v^2)*cos(Theta)         ,  ( c*(u^2+v^2) - w*(a*u+b*v) ) * (1-cos(Theta)) + (a*v-b*u)*sin(Theta)
            0                                  ,  0                                  ,   0                                 ,  1                                                                     ];

ANS = rotMat * [ PointPreRot ; 1 ];

PointAfterRot = ANS(1:3);

end