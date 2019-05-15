%%%%%%%%%%%%%%%%%%%%%%  Index_Local_To_Global_3D  %%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2016-12-04
% Modified on: 2016-12-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function I = Index_Local_To_Global_3D(NX, NY, NZ, i, j, k)

if (i>NX),   error('Current "i" index exceeds NX.');  end
if (j>NY),   error('Current "j" index exceeds NY.');  end
if (k>NZ),   error('Current "k" index exceeds NZ.');  end

I = (k-1)*NX*NY + (j-1)*NX + i;

end