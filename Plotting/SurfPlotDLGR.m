%Plotting - ADM surf plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=linspace(Lx/(2*Grid.Nx), (2*Grid.Nx*Lx-Lx)/(2*Grid.Nx), Grid.Nx);
y=linspace(Ly/(2*Grid.Ny), (2*Grid.Ny*Ly-Ly)/(2*Grid.Ny), Grid.Ny);

% get the corners of the domain in which the data occurs.
min_x = min(min(X));
min_y = min(min(Y));
max_x = max(max(X));
max_y = max(max(Y));
 
% the image data you want to show as a plane.
Nf = FineGrid.Nx*FineGrid.Ny;
for i=1:Nf
    if (CoarseGrid(1).Active(FineGrid.Father(i)) == 1)
        ActiveFine(FineGrid.I(i), FineGrid.J(i)) = 5;
    end
end
planeimg = ActiveFine;
 
% scale image between [0, 255] in order to use a custom color map for it.
minplaneimg = min(min(planeimg)); % find the minimum
scaledimg = (floor(((planeimg - minplaneimg) ./ ...
    (max(max(planeimg)) - minplaneimg)) * 155)); % perform scaling
 
% convert the image to a true color image with the jet colormap.
colorimg = ind2rgb(scaledimg, jet);
 
% set hold on so we can show multiple plots / surfs in the figure.
figure; hold on;
 
% do a normal surface plot.
surf(X,Y,P');
title('Pressure [Pa]');
 
% set a colormap for the surface
colormap(jet);
 
% desired z position of the image plane.
imgzposition = -1e5;
 
% plot the image plane using surf.
colorimg(5:5:end,:,:) = 0;       %# Change every 5 row to white
colorimg(:,5:5:end,:) = 0;       %# Change every 5 column to white 
colorimg(15:15:end,:,:) = 1;       %# Change every 15 row to black
colorimg(:,15:15:end,:) = 1;       %# Change every 15 column to black
surf([min_x max_x],[min_y max_y], repmat(imgzposition, [2 2]),...
    colorimg, 'facecolor','texture');

grid on;
set(gca,'layer','top');
set(gca,'fontsize',24);
% set the view.
view(45,30);
 
% label the axes
xlabel('x');
ylabel('y');
zlabel('Pressure [Pa]');

% set hold on so we can show multiple plots / surfs in the figure.
figure; hold on;
 
% do a normal surface plot.
surf(X,Y,S');
title('Water Saturation');

% set a colormap for the surface
colormap(jet);
 
% desired z position of the image plane.
imgzposition = -1;
 
% plot the image plane using surf.
surf([min_x max_x],[min_y max_y], repmat(imgzposition, [2 2]),...
    colorimg, 'facecolor','texture');

grid on;
set(gca,'layer','top');
set(gca,'fontsize',24);

% set the view.
view(45,30);
 
% label the axes
xlabel('x');
ylabel('y');
zlabel('Water Saturation');