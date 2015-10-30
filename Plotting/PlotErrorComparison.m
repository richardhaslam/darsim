%Plotting - Error contour plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA
DS = 'DS02';
Problem = 'DARSim1';
Nx = 216;
Ny = 54;
Lx = 2160;
Ly = 540;
x=linspace(Lx/(2*Nx), (2*Nx*Lx-Lx)/(2*Nx), Nx);
    y=linspace(Ly/(2*Ny), (2*Ny*Ly-Ly)/(2*Ny), Ny);
    [X, Y] = meshgrid(x,y);
Tcomparison = 10; %goes from 1 to 10
% Read data from files
%Pressure
file = fopen(strcat('../ResultsAndImages/PressureConstrained/',Problem, '/FIM/FIMPressure.txt'), 'r');
dummy = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
FineGrid.P = dummy';
fclose(file);
file = fopen(strcat('../ResultsAndImages/PressureConstrained/',Problem, '/ADM/',DS,'/Constant_2Levels/Pressure.txt'), 'r');
dummy = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
Const2L.P = dummy';
fclose(file);
file = fopen(strcat('../ResultsAndImages/PressureConstrained/',Problem, '/ADM/', DS, '/Bilinear_2Levels/Pressure.txt'), 'r');
dummy = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
Bilin2L.P = dummy';
fclose(file);
%Saturation
file = fopen(strcat('../ResultsAndImages/PressureConstrained/',Problem, '/FIM/FIMSaturation.txt'), 'r');
dummy = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
FineGrid.S = dummy';
fclose(file);
file = fopen(strcat('../ResultsAndImages/PressureConstrained/',Problem, '/ADM/', DS, '/Constant_2Levels/Saturation.txt'), 'r');
dummy = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
Const2L.S = dummy';
fclose(file);
file = fopen(strcat('../ResultsAndImages/PressureConstrained/',Problem, '/ADM/', DS, '/Bilinear_2Levels/Saturation.txt'), 'r');
dummy = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
Bilin2L.S = dummy';
fclose(file);

%Compute Error
 Const2L.errorP =  abs(Const2L.P -  FineGrid.P)./FineGrid.P;
 Bilin2L.errorP =  abs(Bilin2L.P -  FineGrid.P)./FineGrid.P;
 Const2L.errorS =  (Const2L.S -  FineGrid.S);
 Bilin2L.errorS =  (Bilin2L.S -  FineGrid.S);
 
 %Plot error at the end
 figure(103)
 pcolor(X, Y, reshape(Const2L.errorP(:,Tcomparison), Nx, Ny)');
 title(strcat(DS, ': C2L: Pressure relative error [Pa]'));
 xlabel('x [m]');
 ylabel('y [m]');
 colormap(jet);
 colorbar;
 axis('image');
 set(gca,'fontsize',24);
 figure(203)
 pcolor(X, Y, reshape(Bilin2L.errorP(:,Tcomparison), Nx, Ny)');
 title(strcat(DS, ': B2L: Pressure relative error [Pa]'));
 xlabel('x [m]');
 ylabel('y [m]');
 colormap(jet);
 colorbar;
 axis('image');
 set(gca,'fontsize',24);
 figure(106)
 pcolor(X, Y,reshape(Const2L.errorS(:,Tcomparison), Nx, Ny)');
 title(strcat(DS, ': C2L: Absolute saturation error'));
 xlabel('x [m]');
 ylabel('y [m]');
 colormap(jet);
 colorbar;
 caxis([-0.5, 0.5]);
 axis('image');
 set(gca,'fontsize',24);
 figure(206)
 pcolor(X, Y,reshape(Bilin2L.errorS(:,Tcomparison), Nx, Ny)');
 title(strcat(DS, ': B2L: Absolute saturation error'));
 xlabel('x [m]');
 ylabel('y [m]');
 colormap(jet);
 colorbar;
 caxis([-0.5, 0.5]);
 axis('image');
 set(gca,'fontsize',24);