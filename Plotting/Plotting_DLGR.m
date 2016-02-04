%Plotting - ADM pressure and saturation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot for 2D problems
x=linspace(Grid.Lx/(2*Grid.Nx), (2*Grid.Nx*Grid.Lx-Grid.Lx)/(2*Grid.Nx), Grid.Nx);
y=linspace(Grid.Ly/(2*Grid.Ny), (2*Grid.Ny*Grid.Ly-Grid.Ly)/(2*Grid.Ny), Grid.Ny);
x=linspace(0, Grid.Lx-10, Grid.Nx);
y=linspace(0, Grid.Ly-10, Grid.Ny);
[X, Y] = meshgrid(x,y);

%Pressure plot
figure(600)
%suptitle([num2str(round(t/(3600*24))) ' days']);
%subplot(2,1,1);
h = pcolor(X,Y,P');
set(h, 'EdgeColor', 'none');
title('Pressure [Pa]');
xlabel('x [m]');
ylabel('y [m]');
colorbar;
caxis ([0,1e5]);
colormap(jet);
axis('image')
set(gca,'fontsize',24);

hold on
%Fine Grids
Nf = FineGrid.Nx*FineGrid.Ny;
for f=1:Nf
    if (FineGrid.Active(f) == 1)
        Xmax = FineGrid.dx*(FineGrid.I(f));
        Ymax = FineGrid.dy*(FineGrid.J(f));
        Xmin = Xmax - FineGrid.dx;
        Ymin = Ymax - FineGrid.dy;
        
        line([Xmin Xmax], [Ymin Ymin], 'Color', 'w', 'LineWidth',1);
        line([Xmin Xmax], [Ymax Ymax], 'Color', 'w', 'LineWidth',1);
        line([Xmin Xmin], [Ymin Ymax], 'Color', 'w', 'LineWidth',1);
        line([Xmax Xmax], [Ymin Ymax], 'Color', 'w', 'LineWidth',1);
    end
end
%Coarse Grid
for x=1:maxLevel
    Nc = CoarseGrid(x).Nx*CoarseGrid(x).Ny;
    for c=1:Nc
        if CoarseGrid(x).Active(c) == 1
            Xmin = FineGrid.dx*(CoarseGrid(x).I(c) - floor((CoarseGrid(x).CoarseFactor(1) - 1)/2))-FineGrid.dx;
            Xmax = FineGrid.dx*(CoarseGrid(x).I(c) + ceil((CoarseGrid(x).CoarseFactor(1) - 1)/2));
            Ymin = FineGrid.dy*(CoarseGrid(x).J(c) - floor((CoarseGrid(x).CoarseFactor(2) - 1)/2))-FineGrid.dy;
            Ymax = FineGrid.dy*(CoarseGrid(x).J(c) + ceil((CoarseGrid(x).CoarseFactor(2) - 1)/2));            
            line([Xmin Xmax], [Ymin Ymin], 'Color', 'w', 'LineWidth',1);
            line([Xmin Xmax], [Ymax Ymax], 'Color', 'w', 'LineWidth',1);
            line([Xmin Xmin], [Ymin Ymax], 'Color', 'w', 'LineWidth',1);
            line([Xmax Xmax], [Ymin Ymax], 'Color', 'w', 'LineWidth',1);
        end
    end
end

%Saturation plot
figure(601)
%subplot(2,1,2);
h = pcolor(X,Y,S');
set(h, 'EdgeColor', 'none');
title('Water Saturation');
xlabel('x [m]');
ylabel('y [m]');
colorbar;
axis('image');
caxis ([0,1]);
colormap(jet);
set(gca,'fontsize',24);

hold on
%Fine Grids
Nf = FineGrid.Nx*FineGrid.Ny;
for f=1:Nf
    if (FineGrid.Active(f) == 1)
        Xmax = FineGrid.dx*(FineGrid.I(f));
        Ymax = FineGrid.dy*(FineGrid.J(f));
        Xmin = Xmax - FineGrid.dx;
        Ymin = Ymax - FineGrid.dy;
        
        line([Xmin Xmax], [Ymin Ymin], 'Color', 'w', 'LineWidth',1);
        line([Xmin Xmax], [Ymax Ymax], 'Color', 'w', 'LineWidth',1);
        line([Xmin Xmin], [Ymin Ymax], 'Color', 'w', 'LineWidth',1);
        line([Xmax Xmax], [Ymin Ymax], 'Color', 'w', 'LineWidth',1);
    end
end
%Coarse Grid
for x=1:maxLevel
    Nc = CoarseGrid(x).Nx*CoarseGrid(x).Ny;
    for c=1:Nc
        if CoarseGrid(x).Active(c) == 1
            Xmin = FineGrid.dx*(CoarseGrid(x).I(c) - floor((CoarseGrid(x).CoarseFactor(1) - 1)/2))-FineGrid.dx;
            Xmax = FineGrid.dx*(CoarseGrid(x).I(c) + ceil((CoarseGrid(x).CoarseFactor(1) - 1)/2));
            Ymin = FineGrid.dy*(CoarseGrid(x).J(c) - floor((CoarseGrid(x).CoarseFactor(2) - 1)/2))-FineGrid.dy;
            Ymax = FineGrid.dy*(CoarseGrid(x).J(c) + ceil((CoarseGrid(x).CoarseFactor(2) - 1)/2));            
            line([Xmin Xmax], [Ymin Ymin], 'Color', 'w', 'LineWidth',1);
            line([Xmin Xmax], [Ymax Ymax], 'Color', 'w', 'LineWidth',1);
            line([Xmin Xmin], [Ymin Ymax], 'Color', 'w', 'LineWidth',1);
            line([Xmax Xmax], [Ymin Ymax], 'Color', 'w', 'LineWidth',1);
        end
    end
end
drawnow


