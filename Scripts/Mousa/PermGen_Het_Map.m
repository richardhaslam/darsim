close all; clear all; clc;
Directory = '../Permeability/';
Nx = 45;
Ny = 45;
Nz = 01;

File = strcat(Directory, 'logk64_0');
FileID = fopen(File, 'r');
Matrix = textscan(FileID, '%s');
Matrix = str2double(Matrix{1,1});
Matrix = reshape(Matrix,[64,64,64]);
fclose(FileID);

K = zeros(Nx*Ny*Nz,1);
for k = 1:Nz
    for j = 1:Ny
        Index = Nx*Ny*(k-1) + Nx*(j-1) + (1:Nx)';
        K(Index) = Matrix(1:Nx,j,k);
    end
end
K = exp(K);

File = strcat(Directory, 'Perm_Het.txt');
delete(File);
fid = fopen(File,'a+');
fprintf(fid, '%1.6e\n',Nx);
fprintf(fid, '%1.6e\n',Ny);
fprintf(fid, '%1.6e\n',Nz);
fprintf(fid, '%1.6e\n',K);
fclose('all');