%Plotting - permeability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last modified: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPermeability(Perm, Grid)
if (Grid.Nx == 1 || Grid.Ny==1)
    figure(4)
    x = linspace(0, Grid.Lx-Grid.dx, Grid.Nx);
    subplot(3,1,1);
    plot(x, log10(Perm));
    %title('Permeability');
    axis([0 Grid.Lx min(min(log10(Perm)))-1 max(max(log10(Perm)))+1]);
    %axis([0 Grid.Lx 0 1]);
    xlabel('x [m]');
    ylabel('Log(K) [m^2]');
    set(gca,'fontsize',24);
    %axis('image');
    figure(5);
else
    %Plot permeability
    x = linspace(0, Grid.Lx-Grid.dx, Grid.Nx);
    y = linspace(0, Grid.Ly-Grid.dy, Grid.Ny);
    [X, Y] = meshgrid(x,y);
    %Use log scale
    K(:,:) = Perm(1,:,:);
    %K = log(K); %ln
    K = log10(K); %Log
    
    figure(200)
    h = pcolor(X,Y,K');
    set(h, 'EdgeColor', 'none');
    title('X and Y Permeability in log scale');
    xlabel('x [m]');
    ylabel('y [m]');
    axis('image');
    colorbar;
    %caxis([min(min(K)), max(max(K))]);
    colormap(jet);
    set(gca,'fontsize',22);
    hold on
    drawnow;
end
end