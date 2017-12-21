Lx = 9.0; %[m]
Ly = 9.0; %[m]
Lf = 5.0; %[m]
marg_x = (Lx - Lf)/2;
marg_y = (Ly - Lf)/2;

Nx = 15;
Ny = 15;
Nz = 1;

Km = 250; % Matrix Permeability [mD]
% af = mean([Lx,Ly]) / mean([Nx,Ny]); % Fracture Aperture [m]
af = mean([Lx,Ly]) / mean([225,225]); % Averaged/Reduced Aperture [m]
Kf = af^2/12 * 1e15; % Permeability of Fracture (using b^2/12) [mD]
Kf = ( (round(225/Ny)-1)*Km + Kf ) / (round(225/Ny));

K = Km*ones(Nx*Ny,1);
for j=1:Ny
    for i=1:Nx
        I = (j-1)*Ny + i;
        if j==round(Ny/2)
            if i>round(Nx*marg_x/Lx) && i<Nx-round(Nx*marg_x/Lx)
                K(I) = Kf;
            end
        end
        if i==round(Nx/2)
            if j>round(Ny*marg_y/Ly) && j<Ny-round(Ny*marg_y/Ly)
                K(I) = Kf;
            end
        end
    end
end

Directory = '../Permeability/';
File = strcat(Directory, 'Perm_2chan_cross.txt');
delete(File);
fid = fopen(File,'a+');
fprintf(fid, '%1.6e\n',Nx);
fprintf(fid, '%1.6e\n',Ny);
fprintf(fid, '%1.6e\n',Nz);
fprintf(fid, '%1.6e\n',K);
fclose('all');