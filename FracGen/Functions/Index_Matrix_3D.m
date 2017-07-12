%%%%%%%%%%%%%%%%%%%%%%%%%%%  Index_Matrix_3D  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2016-12-04
% Modified on: 2016-12-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function Im = Index_Matrix_3D(NX, NY, NZ, i, j, k)

if (i>NX),   error('Current "i" index exceeds NX.');  end
if (j>NY),   error('Current "j" index exceeds NY.');  end
if (k>NZ),   error('Current "k" index exceeds NZ.');  end

Im = (k-1)*NX*NY + (j-1)*NX + i;

end