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
Geometry = cartGrid([Nx, Ny, Nz], [Lx, Ly, Lz]);                                  % Create Cartesian Grid
Geometry = computeGeometry(Geometry);                                                    % Compute Geomtery of the Cartesian Grid

% SECTION 1: CELLS NODES (X, Y, Z COORDINATES)
Section1 = [(1:1:Geometry.nodes.num)' Geometry.nodes.coords]; 

% SECTION 2: CELLS NODES + CELLS CENTROIDS + CELL VOLUMES
FtC_1 = [gridCellNo(Geometry) Geometry.cells.faces];                                               % Cell + Face + Tag
LI = (FtC_1(:,3) == 1)|(FtC_1(:,3) == 2)|(FtC_1(:,3) == 3)|(FtC_1(:,3) == 4);                               % Identify N S E W faces (1 2 3 4 Tags)
FtC_1(LI,:) = [];                                                                                  % Delete N S E W faces (1 2 3 4 Tags)                            
NtF_1 = [rldecode(1:Geometry.faces.num, diff(Geometry.faces.nodePos), 2) .' Geometry.faces.nodes]; % Face + Nodes
id_faces = ismember(NtF_1(:,1),FtC_1(:,2));                                                                     % Identify the T/B Faces in [Face + Nodes]
NtF_1 = NtF_1(id_faces,:);                                                                               % T/B Faces with their nodes

F = vec2mat(NtF_1(:,1),4);                                                                       % T/B Faces
N = vec2mat(NtF_1(:,2),4);                                                                        % T/B Nodes
NtF_2 = [F(:,1) N];                                                                              % T/B Faces + Nodes

[X,Y] = ismember(FtC_1,NtF_2(:,1));                                                                   % Obtain Index and Location                                      
FNtF_2 = NtF_2(Y(X),:);                                                                              
NtC = reshape(FNtF_2(:,2:5)',8,[])';                                                               % Reshape: Cells vs Nodes
NtC = [NtC(:,4) NtC(:,3) NtC(:,1:2) NtC(:,8) NtC(:,7) NtC(:,5:6)];

faces = Geometry.cells.faces(:,1);
pos   = Geometry.cells.facePos;
cells = 1:Geometry.cells.num ;
FtC_2 = ObtainFaceIndices(faces, pos, cells);                                                              % c = cells

% Cells to Cell
CtC_1 = [gridCellNo(Geometry) Geometry.faces.neighbors(Geometry.cells.faces(:,1),1) Geometry.faces.neighbors(Geometry.cells.faces(:,1),2)];
LI = CtC_1(:,1)==CtC_1(:,2); CtC_1(LI,2) = 0;
LI = CtC_1(:,1)==CtC_1(:,3); CtC_1(LI,3) = 0;
CtC_1 = [CtC_1(:,1) (CtC_1(:,2)+CtC_1(:,3))];
LI = (CtC_1(:,2) == 0);                                                                                % Delete Cells = 0: External Faces
CtC_1(LI,:) = []; 

subs = CtC_1(:,1);
cpc = [1 ; accumarray(subs,1)];
cpc = cumsum(cpc);
    
faces = CtC_1(:,2);
pos = cpc;
CtC_2 = ObtainCellNeighbors(faces, pos, cells);                                                         % c = cells
CtC_2 = sort(CtC_2,2);

Section2 = [(1:1:Geometry.cells.num)' NtC Geometry.cells.centroids Geometry.cells.volumes FtC_2 CtC_2]; % 8 (nodes number per cell + tag nodes)

% SECTION 3: INTERNAL FACES DATA (FACES CONNECTED TO CELLS)
nodes = Geometry.faces.nodes;
pos   = Geometry.faces.nodePos;
faces = 1:Geometry.faces.num ;
NtF_3 = ObtainNodeIndices(nodes, pos, faces);                             % f = Faces

NF = linspace(1, Geometry.faces.num, Geometry.faces.num)';                                      % Create Face Index Vector (Total Number of Faces)

% Assembly Matrix with Face Geometry Data
IF = [NF, Geometry.faces.neighbors, Geometry.faces.areas, Geometry.faces.centroids, Geometry.faces.normals NtF_3];
LI = (IF(:,2) == 0)|(IF(:,3) == 0);                                        % Delete External Faces
IF(LI,:) = [];                                                             % Delete External Faces
% Create Centroid Vector: Face Centroid - Cell Centroid
c_vec = [Geometry.faces.centroids(IF(:,1),:) - Geometry.cells.centroids(IF(:,2),:), Geometry.faces.centroids(IF(:,1),:) - Geometry.cells.centroids(IF(:,3),:)];

% Internal Faces Data: Face Index + Face Area + Face Centroid + Face Normal + Cell Neighbor + Centroid Vector
Section3 = [IF(:,1) IF(:,4:10) IF(:,2) c_vec(:,1:3) IF(:,3) c_vec(:,4:6) IF(:,11:end)];

% SECTION 4: EXTERNAL FACES (FACES AT THE EXTERNAL BOUNDARIES OF THE GRID)
% Assembly Matrix with Face Geometry Data
EF = [NF, Geometry.faces.areas, Geometry.faces.centroids, Geometry.faces.normals, Geometry.faces.neighbors NtF_3];
LI = (EF(:,9) ~= 0)&(EF(:,10) ~= 0);                                       % Delete Internal Faces
EF(LI,:) = [];                                                             % Delete Internal Faces
NC = EF(:,9)+ EF(:,10);                                                                              % Delete Cell Neighboors  == 0

% External Faces Data: Face Index + Face Area + Face Centroid + Face Normal + Cell Neighbor + Centroid Vector
Section4 = [EF(:,1:8) NC (Geometry.faces.centroids(EF(:,1),:) - Geometry.cells.centroids(NC)) EF(:,11:end)];

%% OUTOUT FILE 1: GRID GEOMETRY DATA
OutputFileName = 'CornerPointGrid_DARSim_InputData_CartesianGrid.txt';
disp( '******************* Writing the data into output text file *********************' );
disp(['Writing into file ', OutputFileName]);

fid = fopen(strcat(InputDirectory,'/',OutputFileName) , 'w+' );
fprintf(fid, '%% **************************************************************************\n');
fprintf(fid, '%% CARTESIAN GRID (UNSTRUCTURED FORMAT DATA) FOR DARSIM 2 RESERVOIR SIMULATOR\n');
fprintf(fid, '%% **************************************************************************\n');
fprintf(fid, '--------------------------------------------------------\n');
fprintf(fid, 'RESERVOIR_GRID_NX\n');
fprintf(fid, '%d\n', Geometry.cartDims(1));
fprintf(fid, 'RESERVOIR_GRID_NY\n');
fprintf(fid, '%d\n', Geometry.cartDims(2));
fprintf(fid, 'RESERVOIR_GRID_NZ\n');
fprintf(fid, '%d\n', Geometry.cartDims(3));
fprintf(fid, '--------------------------------------------------------\n');
fprintf(fid, 'ACTIVE_CELLS\n');
fprintf(fid, '%d\n', Geometry.cells.num);
fprintf(fid, '--------------------------------------------------------\n');
fprintf(fid, 'N_FACES\n');
fprintf(fid, '%d\n', Geometry.faces.num);
fprintf(fid, 'N_INTERNAL_FACES\n');
fprintf(fid, '%d\n', size(Section3,1));
fprintf(fid, 'N_EXTERNAL_FACES\n');
fprintf(fid, '%d\n', size(EF,1));
fprintf(fid, '--------------------------------------------------------\n');
fprintf(fid, 'N_NODES\n');
fprintf(fid, '%d\n', Geometry.nodes.num);
fprintf(fid, '--------------------------------------------------------\n');

fprintf(fid, '\n');
fprintf(fid, '%% Section 1: Nodes Coordinates\n');
fprintf(fid, '%% [Node Index] + [Nodes Coordinates (x;y;z)]\n');
fprintf(fid, '\n');
fprintf(fid, 'NODES_COORDINATES\n');
fprintf(fid,'%s %4s %7s %7s\n','Node No.  ','x','y','z');
for ii = 1: size(Section1,1)
    fprintf(fid,'%7d,    %0.3f , %0.3f , %0.3f\n', Section1(ii,:)');
end

fprintf(fid, '\n');
    fprintf(fid, '----------------------------------------------------------------------------------------------------------\n');
    fprintf(fid, '%% Section 2: Cell Data\n');
    fprintf(fid, '%% [Cell Index] + [Cell Nodes] + [Cell Centroid(x;y;z)] + [Cell Volume] + [Faces to Cell] + [Cells to Cell]\n');
    fprintf(fid, '\n');
    fprintf(fid, 'CELL_GEOMETRY\n');
    fprintf(fid,'%s %11s %6s %6s %6s %6s %6s %6s %6s %29s %22s %33s %44s\n','Cell No.  ',...
                'SW_B','SE_B','NE_B','NW_B','SW_T','SE_T','NE_T','NW_T','Cell Centroid(x;y;z)','Cell Volume','Face Indices', 'Cell Neighbor Indices');
            
    FormatSpec = "%7d    ,    %6d,%6d,%6d,%6d,%6d,%6d,%6d,%6d    ,    %6.3f,%6.3f,%6.3f    ,    %12.3f    ,    ";
    for n = 1 : size(Section2,2) - 13 - size(CtC_2,2)
        FormatSpec = strcat(FormatSpec,"%6d,");
    end
    FormatSpec = char(FormatSpec);
    FormatSpec(end)=[];
    FormatSpec = strcat(FormatSpec,"    ,    ");
    for n = 1 : size(Section2,2) - 13 - size(FtC_2,2)
        FormatSpec = strcat(FormatSpec,"%6d,");
    end
    FormatSpec = char(FormatSpec);
    FormatSpec(end)=[];
    
    FormatSpec = strcat(FormatSpec,'\n');
    for ii = 1:size(Section2,1)
        fprintf(fid,FormatSpec, Section2(ii,:)');
    end

fprintf(fid, '\n');
fprintf(fid, '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf(fid, '%% Section 3: Internal Faces (shared with two cells)\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell 1] + [Centroid Vector 1(x,y,z)] + [Neighboring Cell 2] + [Centroid Vector 2(x,y,z)] + [Node Indices]\n');
fprintf(fid, '%% NC: Neighbor Cell\n');
fprintf(fid, '\n');
fprintf(fid, 'INTERNAL_FACE_GEOMETRY\n');
fprintf(fid,'%s %18s %39s %47s %19s %33s %18s %33s %33s\n',...
            'Face No.','Face Area','Face Centroid(x;y;z)','Face Normal(x;y;z)','NC1','Centroid Vector1(x:y:z)','NC2','Centroid Vector2(x,y,z)','Node Indices');

FormatSpec = "%7.0d    ,    %12.3f    ,    %12.3f,%12.3f,%12.3f    ,   %12.3f,%12.3f,%12.3f    ,   %6.0d,%12.3f,%12.3f,%12.3f    ,   %6.0d,%12.3f,%12.3f,%12.3f    ,   ";
for n = 1 :size(Section3,2) - 16
    FormatSpec = strcat(FormatSpec,"%6d,");
end
FormatSpec = char(FormatSpec);
FormatSpec(end)=[];
FormatSpec = strcat(FormatSpec,'\n');
for ii = 1:size(Section3,1)
    fprintf(fid,FormatSpec, Section3(ii,:)');
end

fprintf(fid, '\n');
fprintf(fid, '------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf(fid, '%% Section 4: External Faces (at boundaries | no shared with cells)\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell] + [Centroid Vector (x,y,z)] + [Node Indices]\n');
fprintf(fid, '%% NC: Neighbor Cell\n');
fprintf(fid, '\n');
fprintf(fid, 'EXTERNAL_FACE_GEOMETRY\n');
fprintf(fid,'%s %18s %39s %47s %19s %33s %33s\n',...
            'Face No.','Face Area','Face Centroid(x;y;z)','Face Normal(x;y;z)','NC','Centroid Vector(x,y,z)','Node Indices');

FormatSpec = "%7.0d    ,    %12.3f    ,    %12.3f,%12.3f,%12.3f    ,   %12.3f,%12.3f,%12.3f    ,   %6.0d,%12.3f,%12.3f,%12.3f    ,   ";
for n = 1 :size(Section4,2) - 12
    FormatSpec = strcat(FormatSpec,"%6d,");
end
FormatSpec = char(FormatSpec);
FormatSpec(end)=[];
FormatSpec = strcat(FormatSpec,'\n');
for ii = 1:size(Section4,1)
    fprintf(fid,FormatSpec, Section4(ii,:)');
end
fclose(fid);
end
%% --------------------------------------------------------------------------
function [f, present] = ObtainNodeIndices(nodes, pos, faces)
   eIX = pos;
   nn  = double(diff([pos(faces), ...
                      pos(faces + 1)], [], 2));
   fn  = double(nodes(mcolon(eIX(faces), eIX(faces + 1) - 1), 1));

   m   = numel(faces);
   n   = max(nn);
   f   = nan([n, m]);

   present           = false([max(nodes), 1]);
   present(fn)       = true;

   node_num          = zeros([max(nodes), 1]);
   node_num(present) = 1 : sum(double(present));

   off = reshape((0 : m - 1) .* n, [], 1);

   f(mcolon(off + 1, off + nn)) = node_num(fn);

   f = f .';
end
%% --------------------------------------------------------------------------
function [c, present] = ObtainFaceIndices(faces, pos, cells)
   eIX = pos;
   nn  = double(diff([pos(cells), ...
                      pos(cells + 1)], [], 2));
   fn  = double(faces(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));

   m   = numel(cells);
   n   = max(nn);
   c   = nan([n, m]);

   present           = false([max(faces), 1]);
   present(fn)       = true;

   node_num          = zeros([max(faces), 1]);
   node_num(present) = 1 : sum(double(present));

   off = reshape((0 : m - 1) .* n, [], 1);

   c(mcolon(off + 1, off + nn)) = node_num(fn);
  
   c = c .';
end
%% --------------------------------------------------------------------------
function [c, present] = ObtainCellNeighbors(faces, pos, cells)
   eIX = pos;
   nn  = double(diff([pos(cells), ...
                      pos(cells + 1)], [], 2));
   fn  = double(faces(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));

   m   = numel(cells);
   n   = max(nn);
   c   = nan([n, m]);

   present           = false([max(faces), 1]);
   present(fn)       = true;

   node_num          = zeros([max(faces), 1]);
   node_num(present) = 1 : sum(double(present));

   off = reshape((0 : m - 1) .* n, [], 1);

   c(mcolon(off + 1, off + nn)) = node_num(fn);
  
   c = c .';
end