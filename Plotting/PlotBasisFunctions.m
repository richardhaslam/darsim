%Plot BF used for a certain timestep
x=linspace(0, Lx-10, Grid.Nx);
y=linspace(0, Ly-10, Grid.Ny);
[X, Y] = meshgrid(x,y);

% Coarse Grid MsP
figure(7)
%suptitle('Multilevel Multiscale basis functions')
%subplot(2,1,1);
pcolor(X, Y, reshape(CoarseGrid(2).MsP(:,3),FineGrid.Nx, FineGrid.Ny)');
colorbar;
colormap('jet');
title('BF of Level 3')
xlabel('x [m]');
ylabel('y [m]');
set(gca,'fontsize',24);
axis('image');

figure(14)
%subplot(2,1,2);
pcolor(X, Y, reshape(CoarseGrid(1).MsP(:,22) + CoarseGrid(1).MsP(:,23) + CoarseGrid(1).MsP(:,24) ...
    + CoarseGrid(1).MsP(:,37) + CoarseGrid(1).MsP(:,38) + CoarseGrid(1).MsP(:,39)...
    + CoarseGrid(1).MsP(:,52) + CoarseGrid(1).MsP(:,53) + CoarseGrid(1).MsP(:,54),FineGrid.Nx, FineGrid.Ny)');
title('BF of Level 1')
xlabel('x [m]');
ylabel('y [m]');
colorbar;
colormap('jet');
axis('image');
set(gca,'fontsize',24);

%Plot DLGR Grid and some BF actually used
figure(16)
%subplot(2,1,1)
h = pcolor(X,Y, P');
set(h, 'EdgeColor', 'none');
colorbar;
colormap('jet');
caxis ([0,1e5]);
title('Pressure [Pa]')
axis('image');
xlabel('x [m]');
ylabel('y [m]');
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
figure(8)
%subplot(2,1,2)
%disp(ADMProl(:,246));
h = pcolor(X, Y, reshape(ADMProl(:,50) + ADMProl(:,189) + ADMProl(:,190) + ADMProl(:,191) + ADMProl(:,195) + ADMProl(:,196) + ADMProl(:,197) + ADMProl(:,244) , FineGrid.Nx, FineGrid.Ny)');
%h = pcolor(X, Y, reshape(ADMProl(:,935) + ADMProl(:,2627), FineGrid.Nx, FineGrid.Ny)');
set(h, 'EdgeColor', 'none');
title('Basis Functions')
axis('image');
xlabel('x [m]');
ylabel('y [m]');
set(gca,'fontsize',24);
colormap('jet');
colorbar;
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