%% EggModel permeability

%K = load('../Permeability/EggModel/PERM1_ECL.INC');
for i=1:100
    fileID = fopen(strcat('../Permeability/EggModel/PERM', num2str(i),'_ECL.INC'), 'r');
    matrix = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);
    temp = strfind(matrix{1}, '/'); % Search a specific string and find all rows containing matches
    EndOfPermx = find(~cellfun('isempty', temp));
    rows = EndOfPermx(1) - 3;
    fileID = fopen(strcat('../Permeability/EggModel/EggPerm', num2str(i),'.txt'), 'w');
    size = [60, 60, 7];
    fprintf(fileID, '%d\n', size);
    Perm = matrix{1};
    for row = 3:rows
        temp = str2double(strsplit(char(Perm(row))));
        if ~isnan(temp(1))
            fprintf(fileID, '%3.5e\n', temp(1:end-1));
        end
    end
    fclose(fileID);
end


