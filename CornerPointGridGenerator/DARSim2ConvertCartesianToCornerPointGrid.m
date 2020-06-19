%% ConvertCartesianToCornerPointGrid for DARSim
% This script receives the geometry input from a Cartesian domain and
% creates CornerPointGrid input data file (txt) for DARSim simulator.
% 
% Author: Mousa HosseiniMehr
% Modified on: 2020/04/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DARSim2ConvertCartesianToCornerPointGrid(InputDirectory, InputFileName)
%% Reading the DARSim Main Input File
InputFile = strcat(InputDirectory,'/',InputFileName);
fileID = fopen(InputFile, 'r');
% Read lines from input file
InputMatrix = textscan(fileID, '%s', 'Delimiter', '\n');
InputMatrix = InputMatrix{1};
fclose(fileID);
% Remove lines which are commented (contain --)
Commented = startsWith(InputMatrix, '--');
InputMatrix(Commented) = {'--'}; % removing the string if it is commented.

% Reading the size of the reservoir
temp = strfind(InputMatrix, 'DIMENS');
index = find(~cellfun('isempty', temp));
if isempty(index)
    error('The keyword "DIMENS" is missing. Please check the input file!\n');
end
Lx = str2double(InputMatrix{index+1});
Ly = str2double(InputMatrix{index+2});
Lz = str2double(InputMatrix{index+3});

% Reading number of grid cells
temp = strfind(InputMatrix, 'SPECGRID');
index = find(~cellfun('isempty', temp));
if isempty(index)
    error('The keyword "SPECGRID" is missing. Please check the input file!\n');
end
Nx = str2double(InputMatrix{index+1});
Ny = str2double(InputMatrix{index+2});
Nz = str2double(InputMatrix{index+3});

%% CARTESIAN GRID GEOMETRY DATA: INPUT FILE GENERATION
% SECTION 1: CELL DATA: NODES + CENTROIDS + VOLUMES
G = cartGrid([Nx, Ny, Nz], [Lx, Ly, Lz]);                                  % Create Cartesian Grid
G = computeGeometry(G);                                                    % Compute Geomtery of the Cartesian Grid
NC = linspace(1, G.cells.num, G.cells.num)';                               % Create Cell Index Vector (Total Number of Cells)
CF = G.cells.faces;                                                        % Faces linked to the Cells
LI = (CF(:,2) ~= 5)&(CF(:,2) ~= 6);                                        % Just consider Top and Bttm Faces of the Cell
CF(LI,:) = [];                                                             % Just consider Top and Bttm Faces of the Cell 
ntf = vec2mat(G.faces.nodes,4);                                            % Reshape the Faces/Nodes Column Vector in a Faces/Nodes Matrix
N = reshape((ntf((CF(:,1)),:))',8,[])';                                    % Reshape the Faces/Nodes Matrix
NC = G.nodes.coords;                                                       % X, Y and Z coordinates of the nodes

% Cell Data: Cell Index + Nodes Coordinates (Top/Bttm (NW NE SW SE)) + Centroids (x, y z coordinates) + Volumes
CellData = [NC NC(N(:,4),:) NC(N(:,3),:) NC(N(:,1),:) NC(N(:,2),:) NC(N(:,8),:) NC(N(:,7),:) NC(N(:,5),:) NC(N(:,6),:) G.cells.centroids G.cells.volumes];

% SECTION 2: INTERNAL FACES DATA (FACES CONNECTED TO CELLS)
NF = linspace(1, G.faces.num, G.faces.num)';                      % Create Face Index Vector (Total Number of Faces)
IF = [NF, G.faces.neighbors, G.faces.areas, G.faces.centroids, G.faces.normals];
LI = (IF(:,2) == 0)|(IF(:,3) == 0);                                        % Delete External Faces
IF(LI,:) = [];                                                             % Delete External Faces
% Create Centroid Vector: Face Centroid - Cell Centroid
c_vec = [G.faces.centroids(IF(:,1),:) - G.cells.centroids(IF(:,2),:), G.faces.centroids(IF(:,1),:) - G.cells.centroids(IF(:,3),:)];
% Internal Faces Data: Face Index + Face Area + Face Centroid + Face Normal + Cell Neighbor + Centroid Vector
IF2 = [IF(:,1) IF(:,4:10) IF(:,2) c_vec(:,1:3) IF(:,3) c_vec(:,4:6)];

% SECTION 3: EXTERNAL FACES DATA (FACES AT THE EXTERNAL BOUNDARIES OF THE GRID)
EF = [NF, G.faces.areas, G.faces.centroids, G.faces.normals, G.faces.neighbors];
LI = (EF(:,9) ~= 0)&(EF(:,10) ~= 0);                                       % Delete Internal Faces
EF(LI,:) = [];                                                             % Delete Internal Faces
EF2 = [EF(:,1:8) (EF(:,9)+EF(:,10))];                                      % Delete Cell Neighboors  == 0
% External Faces Data: Face Index + Face Area + Face Centroid + Face Normal + Cell Neighbor + Centroid Vector
EF3 = [EF2 G.faces.centroids(EF2(:,1),:) - G.cells.centroids(EF2(:,9))];

%% OUTOUT FILE 1: GRID GEOMETRY DATA
OutputFileName = 'CornerPointGrid_DARSim_InputData_CartesianGrid.txt';
disp( '******************* Writing the data into output text file *********************' );
disp(['Writing into file ', OutputFileName]);

fid = fopen(strcat(InputDirectory,'/',OutputFileName) , 'w+' );
fprintf(fid, '%% Node Locations for each one of the cells, Nx * Ny *Nz\n');
fprintf(fid, '%% Coordinates X,Y,Z of the eight nones for each one of the cells\n');
fprintf(fid, '\n');
fprintf(fid, '%% The Grid Resolution of the Reservoir:\n');
fprintf(fid, 'RESERVOIR_GRID_NX\n');
fprintf(fid, '%d\n', Nx);
fprintf(fid, 'RESERVOIR_GRID_NY\n');
fprintf(fid, '%d\n', Ny);
fprintf(fid, 'RESERVOIR_GRID_NZ\n');
fprintf(fid, '%d\n', Nz);
fprintf(fid, 'DIMENS_X\n');
fprintf(fid, '%d\n', Lx);
fprintf(fid, 'DIMENS_Y\n');
fprintf(fid, '%d\n', Ly);
fprintf(fid, 'DIMENS_Z\n');
fprintf(fid, '%d\n', Lz);
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
 
fprintf(fid,'%s %28s %30s %32s %30s %31s %31s %31s %31s %30s %20s\n', ...
          'Cell No.  ','North-West Top Corner(x;y;z)','North-East Top Corner(x;y;z)','South-West Top Corner(x;y;z)',...
          'South-East Top Corner(x;y;z)','North-West Bttm Corner(x;y;z)','North-East Bttm Corner(x;y;z)',...
          'South-West Bttm Corner(x;y;z)','South-East Bttm Corner(x;y;z)','Cell Centroid(x;y;z)','Cell Volume');   
     
for ii = 1:size(CellData,1)
    fprintf(fid,'%6.0d , %9d,%9d,%9d , %9d,%9d,%9d , %9d,%9d,%9d , %9d,%9d,%9d , %9d,%9d,%9d , %9d,%9d,%9d , %9d,%9d,%9d , %9d,%9d,%9d , %9d,%9d,%9d , %15.3f\n', CellData(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 2: Faces Connected to Cells\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell 1] + [Centroid Vector 1(x,y,z)] + [Neighboring Cell 2] + [Centroid Vector 2(x,y,z)]\n');
fprintf(fid, '\n');
fprintf(fid, 'INTERNAL_FACE_GEOMETRY\n');
fprintf(fid,'%s %15s %31s %39s %12s %34s %12s %34s\n','Faces No.','Face Area',...
            'Face Centroid(x;y;z)','Face Normal(x;y;z)','NC1','Centroid Vector1(x:y:z)','NC2','Centroid Vector2(x,y,z)');      

for ii = 1:size(IF2,1)
    fprintf(fid,'%8.0d , %13.3f , %11.3f,%11.3f,%11.3f , % 13.3f,%12.3f,% 8.3f , %6.0d , % 12.3f,% 12.3f,% 11.3f , %6.0d , % 12.3f,% 12.3f,% 11.3f\n', IF2(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 3: External Faces (At Boundaries | No Shared with Cells)\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell] + [Centroid Vector (x,y,z)]:\n');
fprintf(fid, '\n');
fprintf(fid, 'EXTERNAL_FACE_GEOMETRY\n');
fprintf(fid,'%s %15s %31s %37s %13s %34s\n','Faces No.','Face Area','Face Centroid(x;y;z)','Face Normal(x;y;z)','NC','Centroid Vector(x,y,z)');

for ii = 1:size(EF3,1)
    fprintf(fid,'%8.0d , %13.3f , %11.3f,%11.3f,%11.3f , % 13.3f,%12.3f,% 8.3f , %6.0d , % 12.3f,% 12.3f, % 11.3f\n', EF3(ii,:)');
end
fclose(fid);

end