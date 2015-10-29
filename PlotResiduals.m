function PlotResiduals(Ro, Rw, Grid)
Nx = Grid.Nx;
Ny = Grid.Ny;
Rp = zeros(Nx*Ny,1);

Rp(:) = Ro(:)+Rw(:);
ro = abs(Ro)/max(abs(Ro));
rw = abs(Rw)/max(abs(Rw));
rp = abs(Rp)/max(abs(Rp));

ro = reshape(ro, Nx, Ny);
rw = reshape(rw, Nx, Ny);
rp=reshape(rp, Nx, Ny);

%Plot Pressure eq. residual
figure(4)
subplot(2,2,1);
pcolor(rp');
if Nx==Ny
    axis square;
end
title(['Residual of pressure eq. ' num2str(max(abs(Rp)))]);
xlabel('x [m]');
ylabel('y [m]');
axis('image');
colorbar;

%Plot oil residual
subplot(2,2,2);
pcolor(ro');
if Nx==Ny
    axis square;
end
title(['Residual of oil eq.: ' num2str(max(abs(Ro)))]);
xlabel('x [m]');
ylabel('y [m]');
axis('image');
colorbar;

%Plot water residual
subplot(2,2,3);
pcolor(rw');
if Nx==Ny
    axis square;
end
title(['Residual of water eq.: ' num2str(norm(Rw,inf))]);
xlabel('x [m]');
ylabel('y [m]');
axis('image');
colorbar;
drawnow;
end