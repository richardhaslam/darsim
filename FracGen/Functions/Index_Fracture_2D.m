%%%%%%%%%%%%%%%%%%%%%%%%%%%  Index_Fracture_2D  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2016-12-04
% Modified on: 2016-12-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function If = Index_Fracture_2D(NX, NY, i, j)

if (i>NX),   error('Current "i" index exceeds NX.');  end
if (j>NY),   error('Current "j" index exceeds NY.');  end

If = (j-1)*NX + i;

end