  %% CornerPointGrid Converter for DARSim2
% Based on MRST (Matlab Reservoir Simulation Toolbox), it reads the input file(in an Eclipse file format)
% and generates the input file with the cell geometry that can be read by DARSim2.
% 
% Author: Janio Piguave
% Modified on: 2020/01/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CORNER POINT GRID GEOMETRY DATA: INPUT FILE GENERATION
close all; clear; clc;
% Read subset of ECLIPSE GRID file. Return: CartDims + ZCORN + COORD + ACTNUM
grdecl = readGRDECL('NPD5.grdecl');
% grdecl = simpleGrdecl([4, 2, 3], 0.12, 'flat', true);
% Extract Corner Point Nodes of the cells (8 nodes * X,Y,Z(coordinates))
[x,y,z] = buildCornerPtNodes(grdecl);
% Dimensions of the Grid in the X, Y and Z coordinates: Vector with the total number of cells
NumberCells = linspace(1, grdecl.cartDims(1)*grdecl.cartDims(2)*grdecl.cartDims(3), grdecl.cartDims(1)*grdecl.cartDims(2)*grdecl.cartDims(3))';
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
C1=[]; c=1;
for k=2:2:size(Z,3)
    Cc(:,:,c)=[Z(:,:,k-1),Z(:,:,k)];                             
    c=c+1;
end
for k=1:size(Z,3)/2
    C1=[C1;Cc(:,:,k)];                                               
end
C1=C1';
C1=reshape(C1,8,[])';

% SECTION 1: CELLS NODES (X, Y, Z) + CENTROIDS + VECTORS A, B, C (X, Y, Z COORDINATES)
% Only Active Cells (Based on ACTNUM info)
Cell_Nodes = [NumberCells, grdecl.ACTNUM, A(:,2), B(:,2), C1(:,2), A(:,6), B(:,6), C1(:,6), A(:,1), B(:,1), C1(:,1),...
              A(:,5), B(:,5), C1(:,5),A(:,4), B(:,4), C1(:,4), A(:,8), B(:,8), C1(:,8), A(:,3), B(:,3), C1(:,3),...
              A(:,7), B(:,7), C1(:,7)];
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
IF2 = [IF(:,1) IF(:,4:10) IF(:,2) c_vec_IF(:,1:3) IF(:,3) c_vec_IF(:,4:6)];

% SECTION 3: EXTERNAL FACES (FACES AT THE EXTERNAL BOUNDARIES OF THE GRID)
EF = [NumberFaces, G.faces.areas, G.faces.centroids, G.faces.normals, G.faces.neighbors];
LI = (EF(:,9) ~= 0)&(EF(:,10) ~= 0);
EF(LI,:) = [];
EF2 = [EF(:,1:8) (EF(:,9)+EF(:,10))];
c_vec_EF = G.faces.centroids(EF2(:,1),:) - G.cells.centroids(EF2(:,9));
EF3 = [EF2 c_vec_EF ];

%% CORNER POINT GRID ROCK PROPERTIES DATA: INPUT FILE GENERATION
ActiveCells = double([NumberCells grdecl.ACTNUM]);
LI = ActiveCells(:,2) == 0;                            % Logical Index                 
ActiveCells(LI,:) = [];                                % Delete Cells that are not active
ActiveCells(:,2) = [];                                 % Delete Columns of Active Cells
% Load Inputs Files
p = load('NPD5_Porosity.txt')';
K = load('NPD5_Permeability.txt')';
% Just Considered Active Cells
p = p(G.cells.indexMap);
poro = p;
poro_text = [ActiveCells poro];
%poro_text2 = [G.cells.indexMap poro];
K = K(G.cells.indexMap);
% Convert K values a diferent units
K = K .* milli * darcy;
perm = bsxfun(@times, [1 1 0.1], K);
perm_txt = [ActiveCells perm];
rock = makeRock(G, perm, poro);
[K, i, j] = permTensor(rock, G.griddim);
K_text = [ActiveCells K];

% %% CALCULATE TRANSMISSIBILITIES
% % Load Inputs Files
% p = load('NPD5_Porosity.txt')';
% K = load('NPD5_Permeability.txt')';
% % Just Considered Active Cells
% p = p(G.cells.indexMap);
% poro = p;
% K = K(G.cells.indexMap);
% % Convert K values a diferent units
% K = K .* milli * darcy;
% perm = bsxfun(@times, [1 1 0.1], K);
% rock = makeRock(G, perm, poro);
% [K, i, j] = permTensor(rock, G.griddim);
% 
% hT1 = zeros(size(IF(:,2)));
% hT2 = zeros(size(IF(:,3)));
% N = IF2(:,6:8);
% C1 = IF2(:,10:12);
% C2 = IF2(:,14:16);
% 
% for k = 1 : size(i,2) 
%     hT1 = hT1 + C1(:, i(k)) .* K(IF(:,2), k) .* N(:, j(k));
% end
% hT1 = hT1 ./ sum(C1.*C1, 2);
% 
% for k = 1 : size(i,2) 
%     hT2 = hT2 + C2(:, i(k)) .* K(IF(:,3), k) .* N(:, j(k));
% end
% hT2 = hT2 ./ sum(C2.*C2, 2);
% hT3 = abs(hT2);
% hTF = 1./(1./hT1 + 1./hT3);
% 
% %Another Approach
% hT_1 = sum(C1(:,:).* N(:,:).* perm(IF(:,2),:), 2) ./ sum(C1(:,:).* C1(:,:),2) ;
% hT_2 = abs(sum(C2(:,:).* N(:,:).* perm(IF(:,3),:), 2) ./ sum(C2(:,:).* C2(:,:),2));
% hT_F = 1./(1./hT_1 + 1./hT_2);
% 
% hT_11 = sum(C1.* N.* perm(IF(:,2),:), 2) ./ sum(C1.* C1,2) ;
% hT_22 = abs(sum(C2.* N.* perm(IF(:,3),:), 2) ./ sum(C2.* C2,2));
% hT_FF = 1./(1./hT_11 + 1./hT_22);
% 
% % hf = IF(:,1);
% % hf2cn1 = IF(:,2);
% % sgn1 = 2*(hf2cn1 == G.faces.neighbors(hf, 1)) - 1;
% % hf2cn2 = IF(:,3);
% % sgn2 = 2*(hf2cn2 == G.faces.neighbors(hf, 1)) - 1;
% % 
% % N1 = IF(:,4:10);
% % N2 = bsxfun(@times, sgn2, IF(:,4:10));

%% OUTOUT FILE: GRID GEOMETRY DATA
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
fprintf(fid, '%d\n', size(IF2,1));
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
      
for ii = 1:size(Cell_Data,1)
    fprintf(fid,'%6.0d , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f;%6.5f;%6.5f , %6.5f\n', Cell_Data(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 2: Faces Connected to Cells\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell 1] + [Centroid Vector 1(x,y,z)] + [Neighboring Cell 2] + [Centroid Vector 2(x,y,z)]\n');
fprintf(fid, '\n');
fprintf(fid, 'INTERNAL_FACE_GEOMETRY\n');
fprintf(fid,'%s %12s %34s %39s %13s %32s %10s %32s\n','Faces No.','Face Area',...
            'Face Centroid(x;y;z)','Face Normal(x;y;z)','NC1','Centroid Vector1(x:y:z)','NC2','Centroid Vector2(x,y,z)');      

for ii = 1:size(IF2,1)
    fprintf(fid,'%8.0d , %13.6f , %11.6f;%11.6f;%11.6f , % 13.6f;%12.6f;% 8.6f , %6.0d , % 12.6f;% 12.6f;% 11.6f , %6.0d , % 12.6f;% 12.6f;% 11.6f  \n', IF2(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 3: External Faces (At Boundaries | No Shared with Cells)\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell] + [Centroid Vector (x,y,z)]:\n');
fprintf(fid, '\n');
fprintf(fid, 'EXTERNAL_FACE_GEOMETRY\n');
fprintf(fid,'%s %12s %34s %39s %13s %32s\n','Faces No.','Face Area','Face Centroid(x;y;z)','Face Normal(x;y;z)','NC','Centroid Vector(x,y,z)');      

for ii = 1:size(EF3,1)
    fprintf(fid,'%8.0d , %13.6f , %11.6f;%11.6f;%11.6f , % 13.6f;%12.6f;% 8.6f , %6.0d , % 12.6f;% 12.6f; % 11.6f \n', EF3(ii,:)');
end
fclose(fid);

%% OUTPUT FILE: ROCK PROPERTIES DATA
Directory = 'C:\Users\Janio Paul\DARSim2\MSRT\mrst-2019a_zip\'; 
OutputFileName = 'CornerPointGrid_DARSim_RockPropertiesData.txt';
disp( '******************* Writing the data into output text file *********************' );
disp(['Writing into file ', OutputFileName]);

fid = fopen(strcat(Directory, OutputFileName), 'w+');
fprintf(fid, '%% Rock Properties of the Johansen Formation\n');
fprintf(fid, '%% Rock Properties Values of Active Cells of the Corner Point Grid Model\n');
fprintf(fid, '%% Permeability tensor is assumed to be diagonal\n');
fprintf(fid, '%% Vertical permeability (Kz) equaling one-tenth of the horizontal permeability (Kx, Ky)\n');
fprintf(fid, '%% Only the x-component Kx is given in the data file. Kx = Ky. Kz = 0.1 Kx\n');
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

fprintf(fid, '\n');
fprintf(fid, '%% Section 1: Porosity\n');
fprintf(fid, '%% [Cell Index] + [Porosity]\n');
fprintf(fid, '%% Units: Fraction\n');
fprintf(fid, 'POROSITY_DATA\n');
fprintf(fid,'%s %s\n', 'Cell No.  ','Porosity');   
      
for ii = 1:size(poro_text,1)
    fprintf(fid,'%6.0d   ,  %.4f\n', poro_text(ii,:)');
end

fprintf(fid, '\n');
fprintf(fid, '%% Section 2: Permeability\n');
fprintf(fid, '%% [Cell Index] + [Permeability (Kx, Ky, Kz)]\n');
fprintf(fid, '%% Permeability Unit (m2 or D or mD)\n');
fprintf(fid, 'PERMEABILITY_UNIT\n');
fprintf(fid, 'm2\n\n');
fprintf(fid, '%% Permeability Scale (Linear or Logarithmic)\n');
fprintf(fid, 'PERMEABILITY_SCALE\n');
fprintf(fid, 'Linear\n\n');
fprintf(fid, 'PERMEABILITY_DATA\n');
fprintf(fid,'%s %7s %10s %10s\n', 'Cell No.  ','Kx','Ky','Kz');   
   
for ii = 1:size(perm_txt,1)
    fprintf(fid,'%6.0d   ,  %.4e;%.4e;%.4e\n', perm_txt(ii,:)');
end

fprintf(fid, '\n');
fprintf(fid, 'PERMEABILITY_TENSOR\n');
 
fprintf(fid,'%s %8s %10s %10s %14s %10s %10s %14s %10s %10s\n', 'Cell No.  ','Kxx','Kxy','Kxz','Kyx','Kyy','Kyz','Kzx','Kzy','Kzz');   
   
for ii = 1:size(K_text,1)
    fprintf(fid,'%6.0d   ,  %.4e;%.4e;%.4e  ,  %.4e;%.4e;%.4e  ,  %.4e;%.4e;%.4e\n', K_text(ii,:)');
end

fclose(fid);