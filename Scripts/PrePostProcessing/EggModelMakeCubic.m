realization = 1;
File = load(strcat('../Permeability/EggModel/EggPerm_ref', num2str(realization),'.txt'), 'w');
K = File(4:end);
N_layer = 60*60;
size = [60, 60, 63];
fileID = fopen(strcat('../Permeability/EggModel/EggPerm_cubic', num2str(realization),'.txt'), 'w');
fprintf(fileID, '%d\n', size);
for i=1:21
    for j=1:3
        z = (i-1)*3 + j;
        K_ref((z-1)*N_layer+1:z*N_layer) = K((i-1)*N_layer+1:i*N_layer);
    end
end
fprintf(fileID, '%3.5e\n', K_ref);
fclose(fileID);