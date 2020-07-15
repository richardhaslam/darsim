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
CF = [gridCellNo(Geometry) Geometry.cells.faces];                                               % Cell + Face + Tag
CF2 = CF;
LI = (CF(:,3) == 1)|(CF(:,3) == 2)|(CF(:,3) == 3)|(CF(:,3) == 4);                               % Identify N S E W faces (1 2 3 4 Tags)
CF(LI,:) = [];                                                                                  % Delete N S E W faces (1 2 3 4 Tags)                            
FN = [rldecode(1:Geometry.faces.num, diff(Geometry.faces.nodePos), 2) .' Geometry.faces.nodes]; % Face + Nodes
F = CF(:,2);                                                                                    % Faces (T B) for each one of the cells
id_f = ismember(FN(:,1),F);                                                                     % Identify the T/B Faces in [Face + Nodes]
FN2 = FN(id_f,:);                                                                               % T/B Faces with their nodes

F2 = vec2mat(FN2(:,1),4);                                                                       % T/B Faces
N = vec2mat(FN2(:,2),4);                                                                        % T/B Nodes
FN3 = [F2(:,1) N];                                                                              % T/B Faces + Nodes

[X,Y] = ismember(F,FN3(:,1));                                                                   % Obtain Index and Location                                      
FN4 = FN3(Y(X),:);                                                                              
FN5 = reshape(FN4(:,2:5)',8,[])';                                                               % Reshape: Cells vs Nodes
FN5 = [FN5(:,4) FN5(:,3) FN5(:,1:2) FN5(:,8) FN5(:,7) FN5(:,5:6)];
FtC = accumarray(CF2(:,1),CF2(:,2),[],@(x){x});                                                 % Find Faces to Cells
FtC2 = padcat(FtC{:});                                                                          % Convert Cell Array to Matrix
FtC2 = FtC2';                                                                                   % Cells X Faces

nf_min = min(diff(Geometry.cells.facePos));                                                     % min # of faces in cells
nf_max = max(diff(Geometry.cells.facePos));                                                     % max # of faces in cells

Section2 = [(1:1:Geometry.cells.num)' FN5 Geometry.cells.centroids Geometry.cells.volumes FtC2]; % 8 (nodes number per cell + tag nodes)

% SECTION 3: INTERNAL FACES DATA (FACES CONNECTED TO CELLS)
% FN = [rldecode(1:Geometry.faces.num, diff(Geometry.faces.nodePos), 2) .' Geometry.faces.nodes]; % Face + Nodes
% 
% NtF = accumarray(FN(:,1),FN(:,2),[],@(x){x});                                                   % Find Nodes to Faces
% NtF2 = padcat(NtF{:});                                                                          % Convert Cell Array to Matrix
% NtF2 = NtF2';                                                                                   % Faces X Nodes

nodes = Geometry.faces.nodes;
pos   = Geometry.faces.nodePos;
faces = 1:Geometry.faces.num ;
NtF2 = get_face_topo(nodes, pos, faces);                             % f = Faces

NF = linspace(1, Geometry.faces.num, Geometry.faces.num)';                                      % Create Face Index Vector (Total Number of Faces)

nn_min = min(diff(Geometry.faces.nodePos));                                                     % min # of nodes in faces
nn_max = max(diff(Geometry.faces.nodePos));                                                     % max # of nodes in faces

% Assembly Matrix with Face Geometry Data
IF = [NF, Geometry.faces.neighbors, Geometry.faces.areas, Geometry.faces.centroids, Geometry.faces.normals NtF2];
LI = (IF(:,2) == 0)|(IF(:,3) == 0);                                        % Delete External Faces
IF(LI,:) = [];                                                             % Delete External Faces
% Create Centroid Vector: Face Centroid - Cell Centroid
c_vec = [Geometry.faces.centroids(IF(:,1),:) - Geometry.cells.centroids(IF(:,2),:), Geometry.faces.centroids(IF(:,1),:) - Geometry.cells.centroids(IF(:,3),:)];

% Internal Faces Data: Face Index + Face Area + Face Centroid + Face Normal + Cell Neighbor + Centroid Vector
IF2 = [IF(:,1) IF(:,4:10) IF(:,2) c_vec(:,1:3) IF(:,3) c_vec(:,4:6) IF(:,11:end)];

% SECTION 4: EXTERNAL FACES (FACES AT THE EXTERNAL BOUNDARIES OF THE GRID)
% Assembly Matrix with Face Geometry Data
EF = [NF, Geometry.faces.areas, Geometry.faces.centroids, Geometry.faces.normals, Geometry.faces.neighbors NtF2];
LI = (EF(:,9) ~= 0)&(EF(:,10) ~= 0);                                       % Delete Internal Faces
EF(LI,:) = [];                                                             % Delete Internal Faces
EF2 = [EF(:,1:8) (EF(:,9)+EF(:,10))];                                      % Delete Cell Neighboors  == 0

% External Faces Data: Face Index + Face Area + Face Centroid + Face Normal + Cell Neighbor + Centroid Vector
EF3 = [EF2 Geometry.faces.centroids(EF2(:,1),:) - Geometry.cells.centroids(EF2(:,9)) EF(:,11:end)];

%% OUTOUT FILE 1: GRID GEOMETRY DATA
OutputFileName = 'CornerPointGrid_DARSim_InputData_CartesianGrid.txt';
disp( '******************* Writing the data into output text file *********************' );
disp(['Writing into file ', OutputFileName]);

fid = fopen(strcat(InputDirectory,'/',OutputFileName) , 'w+' );
fprintf(fid, '%% Node Locations for each one of the cells, Nx * Ny *Nz\n');
fprintf(fid, '%% Coordinates X,Y,Z of the eight nones for each one of the cells\n');
fprintf(fid, '%% NW_T: Northwest Top Corner\n');
fprintf(fid, '%% NE_T: Northeast Top Corner\n');
fprintf(fid, '%% SW_T: Southwest Top Corner\n');
fprintf(fid, '%% SE_T: Southeast Top Corner\n');
fprintf(fid, '%% NW_B: Northwest Bottom Corner\n');
fprintf(fid, '%% NE_B: Northeast Bottom Corner\n');
fprintf(fid, '%% SW_B: Southwest Bottom Corner\n');
fprintf(fid, '%% SE_B: Southeast Bpttom Corner\n');
fprintf(fid, '\n');
fprintf(fid, '%% The Grid Resolution of the Reservoir:\n');
fprintf(fid, 'RESERVOIR_GRID_NX\n');
fprintf(fid, '%d\n', Geometry.cartDims(1));
fprintf(fid, 'RESERVOIR_GRID_NY\n');
fprintf(fid, '%d\n', Geometry.cartDims(2));
fprintf(fid, 'RESERVOIR_GRID_NZ\n');
fprintf(fid, '%d\n', Geometry.cartDims(3));
fprintf(fid, 'ACTIVE_CELLS\n');
fprintf(fid, '%d\n', Geometry.cells.num);
fprintf(fid, 'N_INTERNAL_FACES\n');
fprintf(fid, '%d\n', size(IF2,1));
fprintf(fid, 'N_EXTERNAL_FACES\n');
fprintf(fid, '%d\n', size(EF,1));
fprintf(fid, 'N_NODES\n');
fprintf(fid, '%d\n', Geometry.nodes.num);

fprintf(fid, '\n');
fprintf(fid, '%% Section 1: Nodes Coordinates\n');
fprintf(fid, 'NODES_COORDINATES\n');
fprintf(fid,'%s %7s %17s %17s\n','Node No.  ','x','y','z');
for ii = 1:size(Section1,1)
    fprintf(fid,'%6d ,   %6f ,   %6f ,   %6f\n', Section1(ii,:)');
end

fprintf(fid, '\n');
fprintf(fid, '%% Section 2: Grid Points Coordinates\n');
fprintf(fid, '%% [Nodes Coordinates (x;y;z)] + [Cell Centroid(x;y;z)] + [Cell Volume (m3)]\n');
fprintf(fid, 'CELL_GEOMETRY\n');
fprintf(fid,'%s %6s %9s %9s %9s %9s %9s %9s %9s %33s %24s %17s\n','Cell No.  ',...
    'NW_B','NE_B','SW_B','SE_B','NW_T','NE_T','SW_T','SE_T','Cell Centroid(x;y;z)','Cell Volume','Faces to Cell');

FormatSpec = "%7d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d,    %6.5f,%6.5f,%6.5f ,    %6.5f,    ";
for n = 1 : size(Section2,2) - 13
    FormatSpec = strcat(FormatSpec,"%d,");
end
FormatSpec = char(FormatSpec);
FormatSpec(end)=[];
FormatSpec = strcat(FormatSpec,'\n');
for ii = 1:size(Section2,1)
    fprintf(fid,FormatSpec, Section2(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 3: Faces Connected to Cells\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell 1] + [Centroid Vector 1(x,y,z)] + [Neighboring Cell 2] + [Centroid Vector 2(x,y,z)]\n');
fprintf(fid, '\n');
fprintf(fid, 'INTERNAL_FACE_GEOMETRY\n');
fprintf(fid,'%s %13s %33s %41s %15s %35s %15s %32s %23s\n','Faces No.','Face Area',...
    'Face Centroid(x;y;z)','Face Normal(x;y;z)','NC1','Centroid Vector1(x:y:z)','NC2','Centroid Vector2(x,y,z)','Nodes to Faces');

FormatSpec = "%8.0d , %13.6f , %11.6f,%11.6f,%11.6f , %13.6f,%12.6f,%8.6f , %6.0d , %12.6f,%12.6f,%11.6f , %6.0d , %12.6f,%12.6f,%11.6f,   ";
for n = 1 : size(IF2,2) - 16
    FormatSpec = strcat(FormatSpec,"%d,");
end
FormatSpec = char(FormatSpec);
FormatSpec(end)=[];
FormatSpec = strcat(FormatSpec,'\n');
for ii = 1:size(IF2,1)
    fprintf(fid,FormatSpec, IF2(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 4: External Faces (At Boundaries | No Shared with Cells)\n');
fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell] + [Centroid Vector (x,y,z)]:\n');
fprintf(fid, '\n');
fprintf(fid, 'EXTERNAL_FACE_GEOMETRY\n');
fprintf(fid,'%s %12s %34s %42s %15s %37s %28s\n','Faces No.','Face Area','Face Centroid(x;y;z)','Face Normal(x;y;z)','NC','Centroid Vector(x,y,z)','Nodes to Faces');

FormatSpec = "%8.0d , %13.6f , %11.6f,%11.6f,%11.6f , %13.6f,%12.6f,%8.6f , %6.0d , %12.6f,%12.6f,%11.6f,   ";
for n = 1 : size(EF3,2) - 12
    FormatSpec = strcat(FormatSpec,"%d,");
end
FormatSpec = char(FormatSpec);
FormatSpec(end)=[];
FormatSpec = strcat(FormatSpec,'\n');
for ii = 1:size(EF3,1)
    fprintf(fid,FormatSpec, EF3(ii,:)');
end
fclose(fid);

end
%% --------------------------------------------------------------------------
function [f, present] = get_face_topo(nodes, pos, faces)
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

   tmp         = isfinite(f);
   nnode       = sum(tmp,1);
   ind         = sub2ind(size(f),nnode,1:size(f,2));
   tmp         = repmat(f(ind),size(f,1),1);
   f(isnan(f)) = tmp(isnan(f));
  
   f = f .';
end