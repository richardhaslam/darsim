%% CornerPointGrid Converter for DARSim2
% Based on MRST (Matlab Reservoir Simulation Toolbox), it reads the input file(in an Eclipse file format)
% and generates the input file with the cell geometry that can be read by DARSim2.
% 
% Author: Janio Piguave
% Modified on: 2020/01/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read subset of ECLIPSE GRID file. Return: CartDims + ZCORN + COORD + ACTNUM
grdecl = readGRDECL('NPD5.grdecl');
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

% Matrix with the Nodes (x, y, z) for the Vectors: SW_B, SE_B, NW_B, SW_T
% A = x, B = y, C = z, 
va = [(A(:,7) - A(:,3)), (B(:,7) - B(:,3)), (C(:,7) - C(:,3))];            % vector a = SE_B - SW_B
vb = [(A(:,4) - A(:,3)), (B(:,4) - B(:,3)), (C(:,4) - C(:,3))];            % vector b = NW_B - SW_B
vc = [(A(:,3) - A(:,1)), (B(:,3) - B(:,1)), (C(:,3) - C(:,1))];            % vector c = NW_B - SW_B

% SECTION 1: CELLS NODES (X, Y, Z) + CENTROIDS + VECTORS A, B, C (X, Y, Z COORDINATES)
% Only Active Cells (Based on ACTNUM info)
Nodes_XYZ = [CellsN, grdecl.ACTNUM, A(:,2), B(:,2), C(:,2), A(:,6), B(:,6), C(:,6), A(:,1), B(:,1), C(:,1),...
    A(:,5), B(:,5), C(:,5),A(:,4), B(:,4), C(:,4), A(:,8), B(:,8), C(:,8), A(:,3), B(:,3), C(:,3),...
    A(:,7), B(:,7), C(:,7), va, vb, vc];
TF = Nodes_XYZ(:,2) == 0;
Nodes_XYZ(TF,:) = []; Nodes_XYZ(:,2) = [];
Nodes_XYZi = [Nodes_XYZ(:,1:25), G.cells.centroids, Nodes_XYZ(:,26:34)];

% SECTION 2: CELL DATA (DX - DY - DZ + CELL VOLUMES)
% dx = SE_B_x - SW_B_x  |  dy = NW_B_y - SW_B_y  |  dz = SW_B_z - SW_T_z
delta = [((A(:,7)-A(:,3))), ((B(:,4)-B(:,3))), ((C(:,3)-C(:,1)))];
S2_1 = [CellsN, grdecl.ACTNUM, delta];
TF = S2_1(:,2) == 0;
S2_1(TF,:) = []; S2_1(:,2) = [];
S2_2 = [S2_1, G.cells.volumes];

% SECTION 3: INTERNAL FACES (FACES CONNECTED TO CELLS)
% Create Face Index Vector (Total Number of Faces)
FacesN = linspace(1, G.faces.num, G.faces.num)';
% Assembly Matrix with all the Faces
FC = [FacesN, G.faces.neighbors, G.faces.areas];
TF = (FC(:,2) == 0)|(FC(:,3) == 0);
FC(TF,:) = [];

% SECTION 4: EXTERNAL FACES (FACES AT THE EXTERNAL BOUNDARIES OF THE GRID)
FD = [FacesN, G.faces.neighbors];
TF = (FD(:,2) ~= 0)&(FD(:,3) ~= 0);
FD(TF,:) = [];
%% OUTPUT FILE: WRITING
Directory = 'C:\Users\Janio Paul\DARSim2\MSRT\mrst-2019a_zip\'; 
OutputFileName = 'Test_3.txt';
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
fprintf(fid, '%% [X,Y,Z Nodes Coordinates (Top & Bottom)] + [Cell Centroids] + [Vectors a,b,c]\n');
fprintf(fid, 'CELL_COORD\n');
fprintf(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', ...
          'Cell No.  ', ... 
          'NW_TOP(x) ', '   NW_TOP(y) ', ' NW_TOP(z) ', '   NE_TOP(x)  ', ' NE_TOP(y) ', '  NE_TOP(z) ',...
          '   SW_TOP(x) ', '  SW_TOP(y) ', '  SW_TOP(z) ', '   SE_TOP(x) ', '  SE_TOP(y) ', '  SE_TOP(z) ', ...
          '   NW_BTM(x) ', '  NW_BTM(y) ', '  NW_BTM(z) ', '   NE_BTM(x) ', '  NE_BTM(y) ', '  NE_BTM(z) ', ...
          '   SW_BTM(x) ', '  SW_BTM(y )', '  SW_BTM(z) ', '   SE_BTM(x) ', '  SE_BTM(y) ', '  SE_BTM(z) ', ...
          ' Cell Centroid(x)', 'Cell Centroid(y)', 'Cell Centroid(z)',...
          ' Vector a(x)', 'Vector a(y)', 'Vector a(z)',...
          ' Vector b(x)', 'Vector b(y)', 'Vector b(z)',...
          ' Vector c(x)', 'Vector c(y)', 'Vector c(z)');      
      
for ii = 1:10%size(Nodes_XYZi,1)
    fprintf(fid,'%6.0d: % 6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f ,   %6.5f;   %6.5f;    %6.5f ,    %6.5f;  %6.5f;   %6.5f  ,  %6.5f;  %6.5f;  %6.5f ,    %6.5f;   %6.5f;   %6.5f\n', Nodes_XYZi(ii,:)');
end

fprintf(fid, '\n\n\n');
fprintf(fid, '%% Section 2: Cells Information\n');
fprintf(fid, '%% [Cell Index] + [dx + dy + dz] + [Cell Volume]:\n');
fprintf(fid, '\n');
fprintf(fid, 'CELL_GEOMETRY\n');
fprintf(fid,'%s %s %s %s %s\n','Cell No.  ', '    dx    ', '     dy      ', ' dz  ', '     Cell Volume     ');      

for ii = 1:10%size(S2_2,1)
    fprintf(fid,'%6.0d:   %f;%f;%f , %f\n', S2_2(ii,:)');
end

fprintf(fid, '\n\n\n');
fprintf(fid, '%% Section 3: Faces Connected to Cells  - Data\n');
fprintf(fid, '%% [Face Index] + [Cell Index] + [Cell Index] + [Cell Area]:\n');
fprintf(fid, '\n');
fprintf(fid, 'INTERNAL_FACES\n');
fprintf(fid,'%s %s %s %8s\n','Faces No.','West/Left','East/Right','Area');      

for ii = 1:10%size(FC,1)
    fprintf(fid,'%8.0d: %5d %9d       %.6f\n', FC(ii,:)');
end

fprintf(fid, '\n\n\n');
fprintf(fid, '%% Section 4: External Faces (No Connected to Cells) - Data\n');
fprintf(fid, '%% [Face Index] + [Cell Index] + [Cell Index]:\n');
fprintf(fid, '\n');
fprintf(fid, 'EXTERNAL_FACES\n');
fprintf(fid,'%s %s %s\n','Faces No.','West/Left','East/Right');      

for ii = 1:10%size(FD,1)
    fprintf(fid,'%8.0d: %5d %9d\n', FD(ii,:)');
end
fclose(fid);