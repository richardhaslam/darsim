%Plotting - permeability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPermeability(Perm, Grid, Lx, Ly)
%Plot permeability
x=linspace(0, Lx-10, Grid.Nx);
y=linspace(0, Ly-10, Grid.Ny);
[X, Y] = meshgrid(x,y);
%Use log scale
K(:,:) = Perm(1,:,:);
K = log(K);

figure(100)
h = pcolor(X,Y,K');
set(h, 'EdgeColor', 'none');
title('X and Y Permeability in log scale');
xlabel('x [m]');
ylabel('y [m]');
axis('image');
colorbar;
colormap(jet);
set(gca,'fontsize',22);
hold on

%Draw Wells
% x_i = Inj.x*Grid.dx - 5; 
% y_i = Inj.y*Grid.dy - 5;
% x_p = Prod.x*Grid.dx - 5; 
% y_p = Prod.y*Grid.dy - 5;
% plot(x_i,y_i, 'o');
% plot(x_p,y_p, 'o');
drawnow;
end