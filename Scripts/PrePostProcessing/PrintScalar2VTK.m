 function PrintScalar2VTK(fileID, scalar, name)
 %Print a scalar in VTK format
 fprintf(fileID, strcat('SCALARS ', name,' float 1\n'));
 fprintf(fileID, 'LOOKUP_TABLE default\n');
 fprintf(fileID,'%d ', scalar);
 end