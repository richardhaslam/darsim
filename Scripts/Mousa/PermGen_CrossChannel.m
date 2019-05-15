% This script produces a permeability map with a cross channel in the
% middle. This permeability field is used as reference to validate
% EDFM/pEDFM results.
%% Input
Lx = 1.00; % Reservoir length in x [m]   
Ly = 1.00; % Reservoir length in y [m]
Lz = 0.01; % Reservoir length in z [m]
Nx = 512;  % Reservoir grid cells in x
Ny = 512;  % Reservoir grid cells in y
Nz = 001;  % Reservoir grid cells in z
Km = 1e-14; % Matrix Permeability [m2]

Lf = 0.5; % Fractures Length [m]
margin_x = (Lx - Lf)/2;
margin_y = (Ly - Lf)/2;
% Permeability of Fractures
% af = mean([Lx,Ly]) / mean([Nx,Ny]); % Averaged/reduced fractures aperture [m]
% Kf = af^2/12; % (using b^2/12) [m2]
Kf = 1e-22;

K = Km*ones(Nx*Ny,1);
for j=1:Ny
    for i=1:Nx
        I = (j-1)*Ny + i;
        if j==round(Ny/2)
            if i>round(Nx*margin_x/Lx) && i<Nx-round(Nx*margin_x/Lx)
                K(I) = Kf;
            end
        end
        if i==round(Nx/2)
            if j>round(Ny*margin_y/Ly) && j<Ny-round(Ny*margin_y/Ly)
                K(I) = Kf;
            end
        end
    end
end

Directory = 'D:\SURFdrive\Simulation\DARSim2\Permeability\';
File = strcat(Directory, 'Perm_crossChannel.txt');
delete(File);
fid = fopen(File,'a+');
fprintf(fid, '%1.6e\n',Nx);
fprintf(fid, '%1.6e\n',Ny);
fprintf(fid, '%1.6e\n',Nz);
fprintf(fid, '%1.6e\n',K);
fclose('all');