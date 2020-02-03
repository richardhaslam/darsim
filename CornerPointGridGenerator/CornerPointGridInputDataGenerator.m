%% CornerPointGrid Converter for DARSim2
% Based on MRST (Matlab Reservoir Simulation Toolbox), it reads the input file(in an Eclipse file format)
% and generates the input file with the cell geometry that can be read by DARSim2.
% 
% Author: Janio Piguave
% Modified on: 2020/01/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;
% Read subset of ECLIPSE GRID file. Return: CartDims + ZCORN + COORD + ACTNUM
grdecl = readGRDECL('NPD5.grdecl');
% grdecl = simpleGrdecl([4, 2, 3], 0.12, 'flat', true);
% Extract Corner Point Nodes of the cells (8 nodes * X,Y,Z(coordinates))
[x,y,z] = buildCornerPtNodes(grdecl);
% Dimensions of the Grid in the X, Y and Z coordinates
Nx = grdecl.cartDims(1); Ny = grdecl.cartDims(2); Nz = grdecl.cartDims(3);
% Vector with the total number of cells
CellsN = linspace(1, Nx*Ny*Nz, Nx*Ny*Nz)';
% Compute grid topology and geometry from pillar grid description
G = processGRDECL(grdecl, 'Verbose', true);
% Compute geometry information (centroids, volumes, areas) of the cells
G = computeGeometry(G);

% Reshape the cell data based on the requirements of the input file for DARSim2: X Coordinates
R = size(x);
X = reshape(permute(reshape(x,R(1),2,[],R(3)),[1,3,2,4]),[],2,R(3));

A=[]; a=1;
for k=2:2:size(X,3)
    Aa(:,:,a)=[X(:,:,k-1),X(:,:,k)];                              
    a=a+1;
end
for k=1:size(X,3)/2
    A=[A;Aa(:,:,k)];                                                 
end
A=A';
A=reshape(A,8,[])';

% Reshape the cell data based on the requirements of the input file for DARSim2: Y Coordinates
S = size(y);
Y = reshape(permute(reshape(y,S(1),2,[],S(3)),[1,3,2,4]),[],2,S(3));
B=[]; b=1;
for k=2:2:size(Y,3)
    Bb(:,:,b)=[Y(:,:,k-1),Y(:,:,k)];                            
    b=b+1;
end
for k=1:size(Y,3)/2
    B=[B;Bb(:,:,k)];                                                
end
B=B';
B=reshape(B,8,[])';

% Reshape the cell data based on the requirements of the input file for DARSim2: Z Coordinates
T = size(z);
Z = reshape(permute(reshape(z,T(1),2,[],T(3)),[1,3,2,4]),[],2,T(3));
C=[]; c=1;
for k=2:2:size(Z,3)
    Cc(:,:,c)=[Z(:,:,k-1),Z(:,:,k)];                             
    c=c+1;
end
for k=1:size(Z,3)/2
    C=[C;Cc(:,:,k)];                                               
end
C=C';
C=reshape(C,8,[])';

% SECTION 1: CELLS NODES (X, Y, Z) + CENTROIDS + VECTORS A, B, C (X, Y, Z COORDINATES)
% Only Active Cells (Based on ACTNUM info)
Nodes_XYZ = [CellsN, grdecl.ACTNUM, A(:,2), B(:,2), C(:,2), A(:,6), B(:,6), C(:,6), A(:,1), B(:,1), C(:,1),...
    A(:,5), B(:,5), C(:,5),A(:,4), B(:,4), C(:,4), A(:,8), B(:,8), C(:,8), A(:,3), B(:,3), C(:,3),...
    A(:,7), B(:,7), C(:,7)];
TF = Nodes_XYZ(:,2) == 0;
Nodes_XYZ(TF,:) = []; Nodes_XYZ(:,2) = [];
Nodes_XYZi = [Nodes_XYZ(:,1:25), G.cells.centroids, G.cells.volumes];

% SECTION 2: INTERNAL FACES (FACES CONNECTED TO CELLS)
% Create Face Index Vector (Total Number of Faces)
FacesN = linspace(1, G.faces.num, G.faces.num)';
% Assembly Matrix with all the Faces
FC = [FacesN, G.faces.neighbors, G.faces.areas, G.faces.centroids, G.faces.normals];
TF = (FC(:,2) == 0)|(FC(:,3) == 0);
FC(TF,:) = [];
vcctfc = [G.faces.centroids(FC(:,1),:) - G.cells.centroids(FC(:,2),:), G.faces.centroids(FC(:,1),:) - G.cells.centroids(FC(:,3),:)];
FC2 = [FC(:,1) FC(:,4:10) FC(:,2) vcctfc(:,1:3) FC(:,3) vcctfc(:,4:6)];

% SECTION 3: EXTERNAL FACES (FACES AT THE EXTERNAL BOUNDARIES OF THE GRID)
FD = [FacesN, G.faces.neighbors, G.faces.areas, G.faces.centroids, G.faces.normals];
TF = (FD(:,2) ~= 0)&(FD(:,3) ~= 0);
FD(TF,:) = [];

%% OUTPUT FILE: WRITING
Directory = 'C:\Users\Janio Paul\DARSim2\MSRT\mrst-2019a_zip\'; 
OutputFileName = 'CornerPointGrid_DARSim_InputData.txt';
disp( '******************* Writing the data into output text file *********************' );
disp(['Writing into file ', OutputFileName]);

fid = fopen(strcat(Directory,OutputFileName) , 'w+' );
fprintf(fid, '%% Node Locations for each one of the cells, Nx * Ny *Nz\n');
fprintf(fid, '%% Coordinates X,Y,Z of the eight nones for each one of the cells\n');
fprintf(fid, '\n');
fprintf(fid, '%% The grid resolution of reservoir: [Nx * Ny *Nz]\n');
fprintf(fid, 'RESERVOIR_GRID   %8.0f x %8.0f x %8.0f\n' , grdecl.cartDims(1) , grdecl.cartDims(2) , grdecl.cartDims(3) );
fprintf(fid, '\n');
fprintf(fid, '%% Section 1: Grid Points Coordinates\n');
fprintf(fid, '%% [X,Y,Z Nodes Coordinates (Top & Bottom)] + [Cell Centroid(x,y,z)] + [Cell Volume]\n');
fprintf(fid, 'CELL_GEOMETRY\n');
fprintf(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', ...
          'Cell No.  ', ... 
          'NW_TOP(x) ', '   NW_TOP(y) ', ' NW_TOP(z) ', '   NE_TOP(x)  ', ' NE_TOP(y) ', '  NE_TOP(z) ',...
          '   SW_TOP(x) ', '  SW_TOP(y) ', '  SW_TOP(z) ', '   SE_TOP(x) ', '  SE_TOP(y) ', '  SE_TOP(z) ', ...
          '   NW_BTM(x) ', '  NW_BTM(y) ', '  NW_BTM(z) ', '   NE_BTM(x) ', '  NE_BTM(y) ', '  NE_BTM(z) ', ...
          '   SW_BTM(x) ', '  SW_BTM(y )', '  SW_BTM(z) ', '   SE_BTM(x) ', '  SE_BTM(y) ', '  SE_BTM(z) ', ...
          ' Cell Centroid(x)', 'Cell Centroid(y)', 'Cell Centroid(z)','  Cell Volume');      
      
for ii = 1:10%size(Nodes_XYZi,1)
    fprintf(fid,'%6.0d: % 6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f ,   %6.5f;   %6.5f;    %6.5f      %6.5f\n', Nodes_XYZi(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 2: Faces Connected to Cells\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Cell Index West/Left] + [Centroid Vector W/L (x,y,z)] + [Cell Index East/Right] + [Centroid Vector E/R (x,y,z)]\n');
fprintf(fid, '\n');
fprintf(fid, 'INTERNAL_FACES\n');
fprintf(fid,'%s %12s %18s %s %s %18s %15s %s %12s %s %s %s %12s %s %s %s\n','Faces No.','Face Area','Face Centroid(x)','Face Centroid(y)','Face Centroid(z)','Face Normal(x)','Face Normal(y)','Face Normal(z)','Cell(W/L)','Centroid Vector(x)','Centroid Vector(y)','Centroid Vector(z)','Cell(E/R)','Centroid Vector(x)','Centroid Vector(y)','Centroid Vector(z)');      

for ii = 1:10%size(FC,1)
    fprintf(fid,'%8.0d: %13.6f   %13.6f  ;  %13.6f  ;  %11.6f   ,   % 13.6f ; % 12.6f  ;  % 8.6f    ,  %5d      % 12.6f  ;  % 14.6f    ;   % 11.6f      ,  %5d      % 12.6f  ;  % 14.6f    ;   % 11.6f  \n', FC2(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 3: External Faces (At Boundaries | No Shared with Cells)\n');
fprintf(fid, '%% [Face Index] + [Cell Index] + [Cell Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)]:\n');
fprintf(fid, '\n');
fprintf(fid, 'EXTERNAL_FACES\n');
fprintf(fid,'%s %s %s %9s %21s %s %s %18s %15s %s\n','Faces No.','Cell(W/L)','Cell(E/R)','Area','Face Centroid(x)','Face Centroid(y)','Face Centroid(z)','Face Normal(x)','Face Normal(y)','Face Normal(z)');      

for ii = 1:10%size(FD,1)
    fprintf(fid,'%8.0d: %5d %10d %16.6f   %13.6f  ;  %13.6f  ;  %11.6f   ,   % 13.6f ; % 12.6f  ;  % 8.6f\n', FD(ii,:)');
end
fclose(fid);