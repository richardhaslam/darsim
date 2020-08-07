%% CornerPointGrid Converter for DARSim2
% Based on MRST (Matlab Reservoir Simulation Toolbox), it reads the input file(in an Eclipse file format)and generates the input file 
% with the cell geometry that can be read by DARSim2.
% Author: Janio Piguave
% Modified on: 2020/01/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DARSim2CornerPointGridInputDataGenerator(Directory,FileName)
%% CORNER POINT GRID GEOMETRY DATA: INPUT FILE GENERATION
% Read subset of ECLIPSE GRID file. Return: CartDims + ZCORN + COORD + ACTNUM
InputFile = strcat(Directory,'\',FileName);
grdecl = readGRDECL(InputFile);
Geometries = processGRDECL(grdecl, 'Verbose', true);                                                       % Compute grid topology and geometry from pillar grid description

for G = 1 : length(Geometries)
    Geometry = computeGeometry(Geometries(G));                                                             % Compute geometry of  cells: centroids, volumes, areas                                                                    
    ACTNUM = grdecl.ACTNUM;                                                                                % Removing the active cells of other geometries (if any) from the raw data (grdecl.ACTNUM)
    for i = [1:G-1,G+1:length(Geometries)]
        ACTNUM(Geometries(i).cells.indexMap)= 0;
    end

    %% SECTION 1.1: CELLS NODES (X, Y, Z COORDINATES)
    Section_1 = [(1:1:Geometry.nodes.num)' Geometry.nodes.coords]; 
    
    %% SECTION 1.2: CELLS NODES + CELLS CENTROIDS + CELL VOLUMES
    FtC_1 = [gridCellNo(Geometry) Geometry.cells.faces];                                                   % Cell + Faces to Cell + Faces Tag
    LI = (FtC_1(:,3) == 1)|(FtC_1(:,3) == 2)|(FtC_1(:,3) == 3)|(FtC_1(:,3) == 4);                          % Logical Index: N/S/E/W faces(Faces Tag: 1 2 3 4)
    FtC_1(LI,:) = [];                                                                                      % Delete N/S/E/W faces (Faces Tag: 1 2 3 4))                            
    NtF_1 = [rldecode(1:Geometry.faces.num, diff(Geometry.faces.nodePos), 2).' Geometry.faces.nodes];      % [Face + Nodes]
    id_faces = ismember(NtF_1(:,1),FtC_1(:,2));                                                            % Identify the T/B Faces in [Face + Nodes]
    NtF_1 = NtF_1(id_faces,:);                                                                             % T/B Faces with their nodes

    F = vec2mat(NtF_1(:,1),4);                                                                             % T/B Faces
    N = vec2mat(NtF_1(:,2),4);                                                                             % T/B Nodes
    NtF_2 = [F(:,1) N];                                                                                    % T/B Faces + Nodes

    [X,Y] = ismember(FtC_1(:,2),NtF_2(:,1));                                                               % Obtain Index and Location                                      
    NtF_2 = NtF_2(Y(X),:);                                                                              
    NtC = reshape(NtF_2(:,2:5)',8,[])';                                                                    % Reshape: Cells vs Nodes
    NtC = [NtC(:,4) NtC(:,3) NtC(:,1:2) NtC(:,8) NtC(:,7) NtC(:,5:6)];

    % Faces to Cells
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
    
    Section_2 = [(1:1:Geometry.cells.num)' NtC Geometry.cells.centroids Geometry.cells.volumes FtC_2 CtC_2];  % 8 (nodes number per cell + tag nodes)
        
    %% SECTION 1.3: INTERNAL FACES (FACES CONNECTED TO CELLS)
    nodes = Geometry.faces.nodes;
    pos   = Geometry.faces.nodePos;
    faces = 1:Geometry.faces.num ;
    NtF_3 = ObtainNodeIndices(nodes, pos, faces);                                                             % f = Faces
    
    NF = linspace(1, Geometry.faces.num, Geometry.faces.num)';                                            % Create Face Index Vector (Total Number of Faces)
        
    % Assembly Matrix with Face Geometry Data
    IF = [NF, Geometry.faces.neighbors, Geometry.faces.areas, Geometry.faces.centroids, Geometry.faces.normals NtF_3];
    LI = (IF(:,2) == 0)|(IF(:,3) == 0);                                                                   % Delete External Faces
    IF(LI,:) = [];                                                                                        % Delete External Faces
    % Create Centroid Vector: Face Centroid - Cell Centroid
    c_vec = [(Geometry.faces.centroids(IF(:,1),:) - Geometry.cells.centroids(IF(:,2),:)) (Geometry.faces.centroids(IF(:,1),:) - Geometry.cells.centroids(IF(:,3),:))];

    % Internal Faces Data: Face Index + Face Area + Face Centroid + Face Normal + Cell Neighbor + Centroid Vector
    Section_3 = [IF(:,1) IF(:,4:10) IF(:,2) c_vec(:,1:3) IF(:,3) c_vec(:,4:6) IF(:,11:end)];

    %% SECTION 1.4: EXTERNAL FACES (FACES AT THE EXTERNAL BOUNDARIES OF THE GRID)
    % Assembly Matrix with Face Geometry Data
    EF = [NF, Geometry.faces.areas, Geometry.faces.centroids, Geometry.faces.normals, Geometry.faces.neighbors NtF_3];
    LI = (EF(:,9) ~= 0)&(EF(:,10) ~= 0);                                                                 % Delete Internal Faces
    EF(LI,:) = [];                                                                                       % Delete Internal Faces
    NC = EF(:,9)+ EF(:,10);                                                                              % Delete Cell Neighboors  == 0

    % External Faces Data: Face Index + Face Area + Face Centroid + Face Normal + Cell Neighbor + Centroid Vector
    Section_4 = [EF(:,1:8) NC (Geometry.faces.centroids(EF(:,1),:) - Geometry.cells.centroids(NC)) EF(:,11:end)];

    %% SECTION 1.5: OUTOUT FILE 1 (GRID GEOMETRY DATA)
    disp( '******************* Writing the data into output text files *********************' );
    OutputFileName = 'CornerPointGrid_DARSim_InputData';
    disp(['Writing into file ', OutputFileName, ' #', num2str(G)]);
    
    fid = fopen(strcat(Directory,'\',OutputFileName,'_',num2str(G),'.txt') , 'w+' );
    fprintf(fid, '%% *******************************************************\n');
    fprintf(fid, '%% CORNER POINT GRID DATA FOR DARSIM 2 RESERVOIR SIMULATOR\n');
    fprintf(fid, '%% *******************************************************\n');
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'RESERVOIR_GRID_NX\n');
    fprintf(fid, '%d\n', Geometry.cartDims(1));
    fprintf(fid, '\n');
    fprintf(fid, 'RESERVOIR_GRID_NY\n');
    fprintf(fid, '%d\n', Geometry.cartDims(2));
    fprintf(fid, '\n');
    fprintf(fid, 'RESERVOIR_GRID_NZ\n');
    fprintf(fid, '%d\n', Geometry.cartDims(3));
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'TOTAL_CELLS\n');
    fprintf(fid, '%d\n', Geometry.cartDims(1)*Geometry.cartDims(2)*Geometry.cartDims(3));
    fprintf(fid, '\n');
    fprintf(fid, 'ACTIVE_CELLS\n');
    fprintf(fid, '%d\n', Geometry.cells.num);
    fprintf(fid, '\n');
    fprintf(fid, 'INACTIVE_CELLS\n');
    fprintf(fid, '%d\n', Geometry.cartDims(1)*Geometry.cartDims(2)*Geometry.cartDims(3) - Geometry.cells.num);
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'N_FACES\n');
    fprintf(fid, '%d\n', Geometry.faces.num);
    fprintf(fid, '\n');
    fprintf(fid, 'N_INTERNAL_FACES\n');
    fprintf(fid, '%d\n', size(Section_3,1));
    fprintf(fid, '\n');
    fprintf(fid, 'N_EXTERNAL_FACES\n');
    fprintf(fid, '%d\n', size(Section_4,1));
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'N_NODES\n');
    fprintf(fid, '%d\n', Geometry.nodes.num);
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'MAX_NUMBER_FACES_TO_CELLS\n');
    fprintf(fid, '%d\n', size(FtC_2,2));
    fprintf(fid, '\n');
    fprintf(fid, 'MAX_NUMBER_CELLS_TO_CELLS\n');
    fprintf(fid, '%d\n', size(CtC_2,2));
    fprintf(fid, '\n');
    fprintf(fid, 'MAX_NUMBER_NODES_TO_FACES\n');
    fprintf(fid, '%d\n', size(NtF_3,2));
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, '\n');
    fprintf(fid, '%% Section 1: Nodes Coordinates\n');
    fprintf(fid, '%% [Node Index] + [Nodes Coordinates (x;y;z)]\n');
    fprintf(fid, '\n');
    fprintf(fid, 'NODES_COORDINATES\n');
    fprintf(fid,'%s %7s %16s %13s\n','Node No.  ','x','y','z');
    for ii = 1: size(Section_1,1)
        fprintf(fid,'%7d,    %6f , %6f , %6f\n', Section_1(ii,:)');
    end
    
    fprintf(fid, '\n');
    fprintf(fid, '----------------------------------------------------------------------------------------------------------\n');
    fprintf(fid, '%% Section 2: Cell Data\n');
    fprintf(fid, '%% [Cell Index] + [Cell Nodes] + [Cell Centroid(x;y;z)] + [Cell Volume] + [Faces to Cell] + [Cells to Cell]\n');
    fprintf(fid, '%% NW_T: Northwest Top Corner\n');
    fprintf(fid, '%% NE_T: Northeast Top Corner\n');
    fprintf(fid, '%% SW_T: Southwest Top Corner\n');
    fprintf(fid, '%% SE_T: Southeast Top Corner\n');
    fprintf(fid, '%% NW_B: Northwest Bottom Corner\n');
    fprintf(fid, '%% NE_B: Northeast Bottom Corner\n');
    fprintf(fid, '%% SW_B: Southwest Bottom Corner\n');
    fprintf(fid, '%% SE_B: Southeast Bpttom Corner\n');
    fprintf(fid, '\n');
    fprintf(fid, 'CELL_GEOMETRY\n');
    fprintf(fid,'%s %11s %6s %6s %6s %6s %6s %6s %6s %36s %29s %25s %154s\n','Cell No.  ',...
                'SW_B','SE_B','NE_B','NW_B','SW_T','SE_T','NE_T','NW_T','Cell Centroid(x;y;z)','Cell Volume','Face Indices', 'Neighbor Indices');
            
    FormatSpec = "%7d    ,    %6d,%6d,%6d,%6d,%6d,%6d,%6d,%6d    ,    %6.5f,%6.5f,%6.5f    ,    %12.5f    ,    ";
    for n = 1 : size(Section_2,2) - 13 - size(CtC_2,2)
        FormatSpec = strcat(FormatSpec,"%6d,");
    end
    FormatSpec = char(FormatSpec);
    FormatSpec(end)=[];
    FormatSpec = strcat(FormatSpec,"    ,    ");
    for n = 1 : size(Section_2,2) - 13 - size(FtC_2,2)
        FormatSpec = strcat(FormatSpec,"%6d,");
    end
    FormatSpec = char(FormatSpec);
    FormatSpec(end)=[];
    
    FormatSpec = strcat(FormatSpec,'\n');
    for ii = 1:size(Section_2,1)
        fprintf(fid,FormatSpec, Section_2(ii,:)');
    end
   
    fprintf(fid, '\n');
    fprintf(fid, '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
    fprintf(fid, '%% Section 3: Internal Faces (shared with two cells)\n');
    fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell 1] + [Centroid Vector 1(x,y,z)] + [Neighboring Cell 2] + [Centroid Vector 2(x,y,z)] + [Node Indices]\n');
    fprintf(fid, '%% NC: Neighbor Cell\n');
    fprintf(fid, '\n');
    fprintf(fid, 'INTERNAL_FACE_GEOMETRY\n');
    fprintf(fid,'%s %18s %39s %47s %22s %33s %18s %33s %27s\n',...
                'Face No.','Face Area','Face Centroid(x;y;z)','Face Normal(x;y;z)','NC1','Centroid Vector1(x:y:z)','NC2','Centroid Vector2(x,y,z)','Node Indices');
    
    FormatSpec = "%7.0d    ,    %12.6f    ,    %12.6f,%12.6f,%12.6f    ,   %12.6f,%12.6f,%12.6f    ,   %6.0d,%12.6f,%12.6f,%12.6f    ,   %6.0d,%12.6f,%12.6f,%12.6f    ,   ";
    for n = 1 :size(Section_3,2) - 16
        FormatSpec = strcat(FormatSpec,"%6d,");
    end
    FormatSpec = char(FormatSpec);
    FormatSpec(end)=[];
    FormatSpec = strcat(FormatSpec,'\n');
    for ii = 1:size(Section_3,1)
        fprintf(fid,FormatSpec, Section_3(ii,:)');
    end
    
    fprintf(fid, '\n');
    fprintf(fid, '------------------------------------------------------------------------------------------------------------------------------------------------\n');
    fprintf(fid, '%% Section 4: External Faces (at boundaries | no shared with cells)\n');
    fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell] + [Centroid Vector (x,y,z)] + [Face Topology]\n');
    fprintf(fid, '%% NC: Neighbor Cell\n');
    fprintf(fid, '\n');
    fprintf(fid, 'EXTERNAL_FACE_GEOMETRY\n');
    fprintf(fid,'%s %18s %39s %47s %22s %33s %29s\n',...
                'Face No.','Face Area','Face Centroid(x;y;z)','Face Normal(x;y;z)','NC','Centroid Vector(x,y,z)','Face Topology');
    
    FormatSpec = "%7.0d    ,    %12.6f    ,    %12.6f,%12.6f,%12.6f    ,   %12.6f,%12.6f,%12.6f    ,   %6.0d,%12.6f,%12.6f,%12.6f    ,   ";
    for n = 1 :size(Section_4,2) - 12
        FormatSpec = strcat(FormatSpec,"%6d,");
    end
    FormatSpec = char(FormatSpec);
    FormatSpec(end)=[];
    FormatSpec = strcat(FormatSpec,'\n');
    for ii = 1:size(Section_4,1)
        fprintf(fid,FormatSpec, Section_4(ii,:)');
    end
    fclose(fid);
    
    %% SECTION 2: VTK PLOTTER FILE GENERATOR
    %% SECTION 2.1: NODE POINTS
    Z = max(Geometry.nodes.coords(:,3)) - Geometry.nodes.coords(:,3);                        % Set Depth from shallow to deep 
    Section_1_VTK = [Geometry.nodes.coords(:,1:2) Z];

    %% SECTION 2.2: CELLS
    NtC = NtC - 1;                                                             % VTK files start to count from 0 to N-1 
    Section_2_VTK = [8 * ones(Geometry.cells.num,1) NtC];                                 % 8 (nodes number per cell + tag nodes)

    %% SECTION 2.3: CELL TYPES
    Section_3_VTK = 11 * ones(Geometry.cells.num,1);                                      % 12 is the VTK reference for Hexahedron
    
    %% OUTOUT FILE: CORNER GRID POINT DATA - VTK
    OutputFileName = 'CornerPointGrid_DARSim_VTK_Plotter';
    disp(['Writing into file ', OutputFileName, ' #', num2str(G)]);
    fid = fopen(strcat(Directory,'\', OutputFileName, '_', num2str(G),'.vtk'), 'w+');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'DARSim2 Reservoir Simulator\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

    fprintf(fid, ['POINTS ' num2str(Geometry.nodes.num) ' double\n']);
    for ii = 1:size(Section_1_VTK,1)
        fprintf(fid,'%f %f %f\n', Section_1_VTK(ii,:)');
    end

    fprintf(fid, '\n');
    fprintf(fid, ['CELLS ' num2str(Geometry.cells.num) ' ' num2str(numel(Section_2_VTK)) '\n']);
    for ii = 1:size(Section_2_VTK,1)
        fprintf(fid,'%d %d %d %d %d %d %d %d %d\n', Section_2_VTK(ii,:)');
    end

    fprintf(fid, '\n');
    fprintf(fid, ['CELL_TYPES ' num2str(Geometry.cells.num) '\n']);
    for ii = 1:size(Section_3_VTK,1)
        fprintf(fid,'%d\n', Section_3_VTK(ii,:)');
    end
    fclose(fid);
    
    %% SECTION 3: CORNER POINT GRID ROCK PROPERTIES DATA: INPUT FILE GENERATION
    %% SECTION 3.1: POROSITY
    poro = grdecl.PORO(Geometry.cells.indexMap);
    poro_text = [Geometry.cells.indexMap, poro];

    %% SECTION 3.2: PERMEABILITY
    if ~strcmp(fieldnames(grdecl),'PERMY')
        grdecl.PERMY = grdecl.PERMX;
    end
    if ~strcmp(fieldnames(grdecl),'PERMZ')
        grdecl.PERMZ = grdecl.PERMY;
    end
    perm = [grdecl.PERMX, grdecl.PERMY, grdecl.PERMZ];
    perm = perm(Geometry.cells.indexMap,:);
    % Convert K values from to Darcy to Metric (if not already in metric)
    perm = perm .* milli * darcy;
    perm_txt = [Geometry.cells.indexMap, perm];
    
    %% SECTION 3.3: OUTPUT FILE 2 (ROCK PROPERTIES DATA)
    OutputFileName = 'CornerPointGrid_DARSim_RockPropertiesData';
    disp(['Writing into file ', OutputFileName, ' #', num2str(G)]);
    
    fid = fopen(strcat(Directory,'\', OutputFileName, '_', num2str(G),'.txt'), 'w+');
    fprintf(fid, '%% Rock Properties of the Norne Field\n');
    fprintf(fid, '%% Rock Properties Values of Active Cells of the Corner Point Grid Model\n');
    fprintf(fid, '%% Permeability tensor is assumed to be diagonal\n');
    fprintf(fid, '%% Vertical permeability (Kz) equaling one-tenth of the horizontal permeability (Kx, Ky)\n');
    fprintf(fid, '%% Only the x-component Kx is given in the data file. Kx = Ky. Kz = 0.1 Kx\n');
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
    
    % Writing the porosity data
    if ~isempty(poro_text)
        fprintf(fid, '\n');
        fprintf(fid, '%% Section 1: Porosity\n');
        fprintf(fid, '%% [Cell Index] + [Porosity]\n');
        fprintf(fid, '%% Units: Fraction\n');
        fprintf(fid, 'POROSITY_DATA\n');
        fprintf(fid,'%s %s\n', 'Cell No.  ','Porosity');
        for ii = 1:size(poro_text,1)
            fprintf(fid,'%6.0d   ,  %.4f\n', poro_text(ii,:)');
        end
    end
    
    % Writing the permeability data
    if ~isempty(perm_txt)
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
    end
    fclose(fid);
    
    %% SECTION 3: Generating the VTK Data
    
end
end
%--------------------------------------------------------------------------
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
   
   % To replace the NaN with dummy values use these 5 lines below
   %tmp         = isfinite(f);
   %nnode       = sum(tmp,1);
   %ind         = sub2ind(size(f),nnode,1:size(f,2));
   %tmp         = repmat(f(ind),size(f,1),1);
   %f(isnan(f)) = tmp(isnan(f));
  
   f = f .';
end

% --------------------------------------------------------------------------
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
   
   % To replace the NaN with dummy values use these 5 lines below
   %tmp         = isfinite(c);
   %nnode       = sum(tmp,1);
   %ind         = sub2ind(size(c),nnode,1:size(c,2));
   %tmp         = repmat(c(ind),size(c,1),1);
   %c(isnan(c)) = tmp(isnan(c));
  
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

   % To replace the NaN with dummy values use these 5 lines below
   %tmp         = isfinite(c);
   %nnode       = sum(tmp,1);
   %ind         = sub2ind(size(c),nnode,1:size(c,2));
   %tmp         = repmat(c(ind),size(c,1),1);
   %c(isnan(c)) = tmp(isnan(c));
  
   c = c .';
end