function DARSim2PlotCornerPointGrid_MATLAB(Directory,FileName)
%% CORNER POINT GRID GEOMETRY DATA: INPUT FILE GENERATION
% Read subset of ECLIPSE GRID file. Return: CartDims + ZCORN + COORD + ACTNUM
InputFile = strcat(Directory,'\',FileName);
grdecl = readGRDECL(InputFile);
Geometries = processGRDECL(grdecl, 'Verbose', true);   

for G = 1 : length(Geometries)
    Geometry = computeGeometry(Geometries(G));                                                             % Compute geometry of  cells: centroids, volumes, areas                                                                    
    ACTNUM = grdecl.ACTNUM;                                                                                % Removing the active cells of other geometries (if any) from the raw data (grdecl.ACTNUM)
    for i = [1:G-1,G+1:length(Geometries)]
        ACTNUM(Geometries(i).cells.indexMap)= 0;
    end

%% 1. PLOT GRIDS
n_f2c = [(1:1:Geometry.cells.num)' diff(Geometry.cells.facePos)];                        % Number of faces 2 cells
f2c_min = min(n_f2c(:,2));                                                 % Minimum number of faces in a grid
f2c_max = max(n_f2c(:,2));                                                 % Minimum number of faces in a grid

% number_faces = Number of Faces in a Cell.
% Change this number between MIN (f2c_min) and MAX (f2c_max)
number_faces = 4;

LI = (n_f2c(:,2) ~= number_faces);                                         % Logical Indexing: Identify all cells with different N number of faces per cell
n_f2c(LI,:) = [];n_f2c(:,2) = [];                                          % Delete Cell with different N number of faces

% "number_faces_cell' contains ALL the Cell Indexing of cells with N number of faces.
% You can plot all the cells or select the cell index of a particular cell(s) for plotting.
number_faces_to_cell = n_f2c;
% number_faces_to_cell = [n_f2c(1) n_f2c(2) n_f2c(3)];

%% 1.1 PLOT ALL THE GRIDS (WITH N NUMBER OF FACES) IN THE SAME GRAPH
figure()
hold on
for i = 1:length(number_faces_to_cell)
    plotGrid(Geometry,number_faces_to_cell(i),'FaceColor',[.8 .8 .8]); view(30,50), axis tight off
end

%% 1.2 VIEW THE GRIDS (WITH N NUMBER OF FACES) IN THE WHOLE GRID
pargs = {'EdgeAlpha'; 0.1; 'EdgeColor'; 'k'};
plotGrid(Geometry,'FaceColor','none', pargs{:});

% %% 1.3 PLOT ALL THE GRIDS WITH A N NUMBER OF FACES IN DIFFERENT GRAPHS
% for i = 1:length(number_faces_to_cell)
%     figure()
%     plotGrid(Geometry,number_faces_to_cell(i),'FaceColor',[.8 .8 .8]); view(30,50), axis tight off
% end

%% 2. Plot Faces
n_n2f = [(1:1:Geometry.faces.num)' diff(Geometry.faces.nodePos)];                        % Number of nodes 2 faces
n2f_min = min(n_n2f(:,2));                                                 % Minimum number of nodes in a face
n2f_max = max(n_n2f(:,2));                                                 % Maximum number of nodes in a face

% number_nodes = Number of Nodes in a Cell.
% Change this number between MIN (n2f_min) and MAX (n2f_max)
number_nodes = 3;
LI = (n_n2f(:,2) ~= number_nodes);                                         % Logical Indexing: Identify all faces with different N number of nodes                                      
n_n2f(LI,:) = []; n_n2f(:,2) = [];                                         % Delete Faces with different N number of nodes

% "number_faces_cell' contains ALL the Cell Indexing of cells with N number of faces.
% You can plot all the cells or select the cell index of a particular cell(s) for plotting.
number_nodes_to_face = n_n2f;
% number_nodes_to_face = [n_n2f(1) n_n2f(2) n_n2f(3)];

%% 2.1 PLOT ALL THE FACES (WITH N NUMBER OF NODES) IN THE SAME GRAPH
figure()
hold on
for i = 1:length(number_nodes_to_face)
    plotFaces(Geometry,number_nodes_to_face(i),'FaceColor',[.8 .8 .8]); view(30,50), axis tight off
end

%% 2.2 VIEW THE GRIDS (WITH N NUMBER OF NODES) IN THE WHOLE GRID
pargs = {'EdgeAlpha'; 0.1; 'EdgeColor'; 'k'};
plotGrid(Geometry,'FaceColor','none', pargs{:});

% %% 2.3 PLOT ALL THE FACES WITH A N NUMBER OF NODES IN DIFFERENT GRAPHS
% for i = 1:length(number_nodes_to_face)
%     figure()
%     plotGrid(Geometry,number_nodes_to_face(i),'FaceColor',[.8 .8 .8]); view(30,50), axis tight off
% end
end