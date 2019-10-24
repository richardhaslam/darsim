% construct permeability

perm = ones(99,2,99);

perm(1:50,:,:) = 1e-12 * ones(50,2,99);
perm(51:99,:,:) = 1e-14 * ones(49,2,99);

fileID = fopen('perm3D.txt','w');
fprintf(fileID,'%6e \n',99);
fprintf(fileID,'%6e \n',2);
fprintf(fileID,'%6e \n',99);
fprintf(fileID,'%6e \n',perm(:));