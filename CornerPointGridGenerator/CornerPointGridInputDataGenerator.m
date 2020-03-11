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
% Reshape the cell data based on the requirements of the input file for DARSim2: X Y Z Coordinates
A = inputdataDARSim(x);
B = inputdataDARSim(y);
C = inputdataDARSim(z);

% SECTION 1: CELLS NODES (X, Y, Z) + CENTROIDS + VECTORS A, B, C (X, Y, Z COORDINATES)
% Only Active Cells (Based on ACTNUM info)
Cell_Nodes = [NumberCells, double(grdecl.ACTNUM), A(:,2), B(:,2), C(:,2), A(:,6), B(:,6), C(:,6), A(:,1), B(:,1), C(:,1),...
              A(:,5), B(:,5), C(:,5),A(:,4), B(:,4), C(:,4), A(:,8), B(:,8), C(:,8), A(:,3), B(:,3), C(:,3),...
              A(:,7), B(:,7), C(:,7)];
% Cell_Nodes = double(Cell_Nodes);
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
EF3 = [EF2 G.faces.centroids(EF2(:,1),:) - G.cells.centroids(EF2(:,9))];

%% Corner Grid Point Data for Plotting - VTK file
% SECTION 1: Point Data for each Active Cell
A_VTK = [double(grdecl.ACTNUM) A];
LI = A_VTK(:,1) == 0;                            % Logical Index                 
A_VTK(LI,:) = [];                                % Delete Cells that are not active
A_VTK(:,1) = [];                                 % Delete Columns of Active Cells
A_VTK = reshape(transpose(A_VTK),1,[])';

B_VTK = [double(grdecl.ACTNUM) B];
LI = B_VTK(:,1) == 0;                            % Logical Index                 
B_VTK(LI,:) = [];                                % Delete Cells that are not active
B_VTK(:,1) = [];                                 % Delete Columns of Active Cells
B_VTK = reshape(transpose(B_VTK),1,[])';

C_VTK = [double(grdecl.ACTNUM) C];
LI = C_VTK(:,1) == 0;                            % Logical Index                 
C_VTK(LI,:) = [];                                % Delete Cells that are not active
C_VTK(:,1) = [];                                 % Delete Columns of Active Cells
C_VTK = reshape(transpose(C_VTK),1,[])';

Cell_VTK1 = [A_VTK B_VTK C_VTK];
% SECTION 2: Number of Points for each Active Cell
Cell_VTK2 = [8 * ones(G.cells.num,1) vec2mat([0 1:(size(Cell_VTK1,1)-1)],8)];

% SECTION 3: Types of Cells: Unstructured = 12
Cell_VTK3 = 11 * ones(G.cells.num,1);
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

%% OUTOUT FILE 1: GRID GEOMETRY DATA
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
    fprintf(fid,'%6.0d , %6.5f,%6.5f,%6.5f , %6.5f,%6.5f,%6.5f , %6.5f,%6.5f,%6.5f , %6.5f,%6.5f,%6.5f , %6.5f,%6.5f,%6.5f , %6.5f,%6.5f,%6.5f , %6.5f,%6.5f,%6.5f , %6.5f,%6.5f,%6.5f , %6.5f,%6.5f,%6.5f , %6.5f\n', Cell_Data(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 2: Faces Connected to Cells\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell 1] + [Centroid Vector 1(x,y,z)] + [Neighboring Cell 2] + [Centroid Vector 2(x,y,z)]\n');
fprintf(fid, '\n');
fprintf(fid, 'INTERNAL_FACE_GEOMETRY\n');
fprintf(fid,'%s %12s %34s %39s %13s %32s %10s %32s\n','Faces No.','Face Area',...
            'Face Centroid(x;y;z)','Face Normal(x;y;z)','NC1','Centroid Vector1(x:y:z)','NC2','Centroid Vector2(x,y,z)');      

for ii = 1:size(IF2,1)
    fprintf(fid,'%8.0d , %13.6f , %11.6f,%11.6f,%11.6f , % 13.6f,%12.6f,% 8.6f , %6.0d , % 12.6f,% 12.6f,% 11.6f , %6.0d , % 12.6f,% 12.6f,% 11.6f\n', IF2(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 3: External Faces (At Boundaries | No Shared with Cells)\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell] + [Centroid Vector (x,y,z)]:\n');
fprintf(fid, '\n');
fprintf(fid, 'EXTERNAL_FACE_GEOMETRY\n');
fprintf(fid,'%s %12s %34s %39s %13s %32s\n','Faces No.','Face Area','Face Centroid(x;y;z)','Face Normal(x;y;z)','NC','Centroid Vector(x,y,z)');

for ii = 1:size(EF3,1)
    fprintf(fid,'%8.0d , %13.6f , %11.6f,%11.6f,%11.6f , % 13.6f,%12.6f,% 8.6f , %6.0d , % 12.6f,% 12.6f, % 11.6f\n', EF3(ii,:)');
end
fclose(fid);

%% OUTOUT FILE 2: CORNER GRID POINT DATA - VTK
Directory = 'C:\Users\Janio Paul\DARSim2\MSRT\mrst-2019a_zip\';
OutputFileName = 'CornerPointGrid_DARSim_VTKData.vtk';
disp( '******************* Writing the data into output text file *********************' );
disp(['Writing into file ', OutputFileName]);
fileID = fopen(strcat(Directory,OutputFileName) , 'w+' );
fprintf(fileID, '# vtk DataFile Version 2.0\n');
fprintf(fileID, 'DARSim 2 Reservoir Simulator\n');
fprintf(fileID, 'ASCII\n');
fprintf(fileID, '\n');
fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');

fprintf(fileID, ['POINTS ' num2str(size(Cell_VTK1,1)) ' double\n']);
for ii = 1:size(Cell_VTK1,1)
    fprintf(fileID,'%f %f %f\n', Cell_VTK1(ii,:)');
end

fprintf(fileID, '\n');
fprintf(fileID, ['CELLS ' num2str(G.cells.num) ' ' num2str(size(Cell_VTK1,1) + G.cells.num) '\n']);
for ii = 1:size(Cell_VTK2,1)
    fprintf(fileID,'%d %d %d %d %d %d %d %d %d\n', Cell_VTK2(ii,:)');
end

fprintf(fileID, '\n');
fprintf(fileID, ['CELL_TYPES ' num2str(G.cells.num) '\n']);
for ii = 1:size(Cell_VTK3,1)
    fprintf(fileID,'%d\n', Cell_VTK3(ii,:)');
end
fclose(fileID);
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
    fprintf(fid,'%6.0d   ,  %.4e,%.4e,%.4e\n', perm_txt(ii,:)');
end

fprintf(fid, '\n');
fprintf(fid, 'PERMEABILITY_TENSOR\n');
 
fprintf(fid,'%s %8s %10s %10s %14s %10s %10s %14s %10s %10s\n', 'Cell No.  ','Kxx','Kxy','Kxz','Kyx','Kyy','Kyz','Kzx','Kzy','Kzz');   

for ii = 1:size(K_text,1)
    fprintf(fid,'%6.0d   ,  %.4e,%.4e,%.4e  ,  %.4e,%.4e,%.4e  ,  %.4e,%.4e,%.4e\n', K_text(ii,:)');
end
fclose(fid);