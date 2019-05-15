%%%%%%%%%%%%%%%%%%%%%%  Index_Global_To_Local_3D  %%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2016-12-04
% Modified on: 2016-12-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [i,j,k] = Index_Global_To_Local_3D(NX, NY, NZ, I)

% Retrieving i,j,k from the "Im" index
i = mod(   I                , NX )   ;   if ( i==0 ),   i = NX;   end
j = mod(  (I-i)/NX          , NY ) +1;   if ( j==0 ),   j = NY;   end
k = mod( ((I-i)/NX -j+1)/NY , NZ ) +1;   if ( k==0 ),   k = NZ;   end
if Index_Local_To_Global_3D(NX, NY, NZ, i, j, k) ~= I
    error('i,j,k are not correspondent with "I" index. Check the formula again!');
end

end