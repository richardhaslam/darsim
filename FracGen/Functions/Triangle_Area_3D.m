%%%%%%%%%%%%%%%%%%%%%%%%  Triangle_Area_3D  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2017-01-02
% Modified on: 2017-01-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = Triangle_Area_3D(P1,P2,P3)

P1P2 = P2 - P1;

P1P3 = P3 - P1;

A = 0.5 * norm( cross ( P1P2 , P1P3 ) );

end