% Average saturation
for fileIndex=1:35
    
    Sol = load(strcat('../Input/AvSatPlot/Output/Solution/AvSatPlot_Sol', num2str(fileIndex-1),'.txt'));
    S = Sol(:, 3);
    
    N = [15, 15];
    cf = [5, 5];
    Nc = N./cf;
    
    R = sparse(prod(Nc), prod(N));
    
    for i=1:Nc(1)
        for j=1:Nc(2)
            indexC = i + Nc(1)*(j-1);
            Imin = (i-1)*cf(1)+1;
            Imax = i*cf(1);
            Jmin = (j-1)*cf(2)+1;
            Jmax = j*cf(2);
            If = Imin:Imax;
            Jf = Jmin:Jmax;
            [p, q] = meshgrid(If, Jf);
            pairs = [p(:), q(:)];
            %indexes of the fine cells
            IndexesF = pairs(:,1) + (pairs(:,2)-1)*N(1);
            R(indexC, IndexesF) = 1;
        end
    end
    S_rest = R * S;
    Sav = R' * (S_rest ./ sum(R, 2));
    
    fileID = fopen(strcat('../Input/AvSatPlot/Output/VTK/AvSatPlot', num2str(fileIndex),'.vtk'), 'a');
    
    %Print a scalar in VTK format
    fprintf(fileID, ' \n');
    name = ' AvSat';
    fprintf(fileID, strcat('SCALARS  ', name,' float 1\n'));
    fprintf(fileID, 'LOOKUP_TABLE default\n');
    %fprintf(fileID,'%d ', scalar);
    fwrite(fileID, Sav','float', 'b');
    
    fclose(fileID);

end