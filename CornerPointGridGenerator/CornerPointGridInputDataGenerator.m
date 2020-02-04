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
NumberCells = linspace(1, Nx*Ny*Nz, Nx*Ny*Nz)';
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
Cell_Nodes = [NumberCells, grdecl.ACTNUM, A(:,2), B(:,2), C(:,2), A(:,6), B(:,6), C(:,6), A(:,1), B(:,1), C(:,1),...
    A(:,5), B(:,5), C(:,5),A(:,4), B(:,4), C(:,4), A(:,8), B(:,8), C(:,8), A(:,3), B(:,3), C(:,3),...
    A(:,7), B(:,7), C(:,7)];
LI = Cell_Nodes(:,2) == 0;                            % Logical Index                 
Cell_Nodes(LI,:) = [];                                % Delete Cells that are not active
Cell_Nodes(:,2) = [];                                 % Delete Columns of Active Cells
Cell_Data = [Cell_Nodes, G.cells.centroids, G.cells.volumes];

% SECTION 2: INTERNAL FACES (FACES CONNECTED TO CELLS)
% Create Face Index Vector (Total Number of Faces)
NumberFaces = linspace(1, G.faces.num, G.faces.num)';
% Assembly Matrix with all the Faces
IF = [NumberFaces, G.faces.neighbors, G.faces.areas, G.faces.centroids, G.faces.normals];
LI = (IF(:,2) == 0)|(IF(:,3) == 0);
IF(LI,:) = [];
c_vec_IF = [G.faces.centroids(IF(:,1),:) - G.cells.centroids(IF(:,2),:), G.faces.centroids(IF(:,1),:) - G.cells.centroids(IF(:,3),:)];
FC2 = [IF(:,1) IF(:,4:10) IF(:,2) c_vec_IF(:,1:3) IF(:,3) c_vec_IF(:,4:6)];

% SECTION 3: EXTERNAL FACES (FACES AT THE EXTERNAL BOUNDARIES OF THE GRID)
EF = [NumberFaces, G.faces.areas, G.faces.centroids, G.faces.normals, G.faces.neighbors];
LI = (EF(:,9) ~= 0)&(EF(:,10) ~= 0);
EF(LI,:) = [];
EF2 = [EF(:,1:8) (EF(:,9)+EF(:,10))];
c_vec_EF = G.faces.centroids(EF2(:,1),:) - G.cells.centroids(EF2(:,9));
EF3 = [EF2 c_vec_EF ];

%% OUTPUT FILE: WRITING
Directory = 'C:\Users\Janio Paul\DARSim2\MSRT\mrst-2019a_zip\'; 
OutputFileName = 'CornerPointGrid_DARSim_InputData.txt';
disp( '******************* Writing the data into output text file *********************' );
disp(['Writing into file ', OutputFileName]);

fid = fopen(strcat(Directory,OutputFileName) , 'w+' );
fprintf(fid, '%% Node Locations for each one of the cells, Nx * Ny *Nz\n');
fprintf(fid, '%% Coordinates X,Y,Z of the eight nones for each one of the cells\n');
fprintf(fid, '\n');
fprintf(fid, '%% The Grid Resolution of the Reservoir:\n');
fprintf(fid, 'RESERVOIR_GRID_NX\n');
fprintf(fid, '%d\n', grdecl.cartDims(1));
fprintf(fid, 'RESERVOIR_GRID_NY\n');
fprintf(fid, '%d\n', grdecl.cartDims(2));
fprintf(fid, 'RESERVOIR_GRID_NZ\n');
fprintf(fid, '%d\n', grdecl.cartDims(3));
fprintf(fid, 'ACTIVE_CELLS\n');
fprintf(fid, '%d\n', G.cells.num);
fprintf(fid, 'N_INTERNAL_FACES\n');
fprintf(fid, '%d\n', size(FC2,1));
fprintf(fid, 'N_EXTERNAL_FACES\n');
fprintf(fid, '%d\n', size(EF,1));

fprintf(fid, '\n');
fprintf(fid, '%% Section 1: Grid Points Coordinates\n');
fprintf(fid, '%% [Nodes Coordinates (x;y;z) (Top & Bottom)] + [Cell Centroid(x;y;z)] + [Cell Volume]\n');
fprintf(fid, 'CELL_GEOMETRY\n');
 
fprintf(fid,'%s %31s %39s %39s %39s %40s %38s %39s %39s %35s %22s\n', ...
          'Cell No.  ','North-West Top Corner(x;y;z)','North-East Top Corner(x;y;z)','South-West Top Corner(x;y;z)',...
          'South-East Top Corner(x;y;z)','North-West Bttm Corner(x;y;z)','North-East Bttm Corner(x;y;z)',...
          'South-West Bttm Corner(x;y;z)','South-East Bttm Corner(x;y;z)','Cell Centroid(x;y;z)','Cell Volume');   
      
for ii = 1:10%size(Nodes_XYZi,1)
    fprintf(fid,'%6.0d , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f\n', Cell_Data(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 2: Faces Connected to Cells\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell 1] + [Centroid Vector 1(x,y,z)] + [Neighboring Cell 2] + [Centroid Vector 2(x,y,z)]\n');
fprintf(fid, '\n');
fprintf(fid, 'INTERNAL_FACE_GEOMETRY\n');
fprintf(fid,'%s %12s %34s %39s %13s %32s %10s %32s\n','Faces No.','Face Area',...
            'Face Centroid(x;y;z)','Face Normal(x;y;z)','NC1','Centroid Vector1(x:y:z)','NC2','Centroid Vector2(x,y,z)');      

for ii = 1:10%size(FC2,1)
    fprintf(fid,'%8.0d , %13.6f , %11.6f;%11.6f;%11.6f , % 13.6f;%12.6f;% 8.6f , %d , % 12.6f;% 12.6f;% 11.6f , %d , % 12.6f;% 12.6f;% 11.6f  \n', FC2(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 3: External Faces (At Boundaries | No Shared with Cells)\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell] + [Centroid Vector (x,y,z)]:\n');
fprintf(fid, '\n');
fprintf(fid, 'EXTERNAL_FACE_GEOMETRY\n');
fprintf(fid,'%s %12s %34s %39s %13s %32s\n','Faces No.','Face Area','Face Centroid(x;y;z)','Face Normal(x;y;z)','NC','Centroid Vector(x,y,z)');      

for ii = 1:10%size(FD,1)
    fprintf(fid,'%8.0d , %13.6f , %11.6f;%11.6f;%11.6f , % 13.6f;%12.6f;% 8.6f , %d , % 12.6f;% 12.6f; % 11.6f \n', EF3(ii,:)');
end
fclose(fid);