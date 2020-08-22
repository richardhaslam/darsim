% Reader for eclipse input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim Reservoir Simulator
% Author: Mousa HosseiniMehr
% TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reader_eclipse < handle
    properties
        Directory
        File
    end
    methods
        function obj = reader_eclipse(dir, file)
            if nargin == 1
                obj.Directory = [];
                obj.File = dir;
            else
                obj.Directory = dir;
                obj.File = strcat(dir,'/',file);
            end
        end
        function CornerPointGridData = ReadGRDECL(obj,Index_List)
            % Read the ECLIPSE GRDECL file
            fprintf('---> Reading GRDECL raw data ... ');
            grdecl_raw = readGRDECL(obj.File);
            fprintf('Done!\n');
            % Process the GRDECL data
            fprintf('---> Pocessing GRDECL raw data ... ');
            grdecl_processed = mprocessGRDECL(grdecl_raw, 'Verbose', true);
            fprintf('Done!\n');
            
            % Some models have more than one disconnected formations. Use "length(Geometries)" as
            % the maximum counter for the following for-loop to process all of them.
            if nargin > 1
                switch Index_List
                    case{'all','All','ALL'}
                        Index_List = 1:length(grdecl_processed);
                    otherwise
                        if any(Index_List<1) || any(Index_List>length(grdecl_processed))
                            error('DARSim Error: The index list for reading the ECLIPSE data is out of bound (should be in range of 1 to %d.',length(grdecl_processed));
                        end
                end
            else
                Index_List = 1;
            end

            % Looping over the list of processed grdecl data (usually there will be only 1)
            CornerPointGridData(1:length(Index_List),1) = struct;
            for n = 1:length(Index_List)
                G = Index_List(n);
                CornerPointGridData(n).Index = G;

                % Remove the active cells of other geometries (if any) from the raw data (grdecl.ACTNUM)
                ACTNUM = grdecl_raw.ACTNUM;
                for i = [1:G-1,G+1:length(grdecl_processed)]
                    ACTNUM(grdecl_processed(i).cells.indexMap)= 0;
                end

                % Compute geometry of cells: centroids, volumes, areas
                fprintf('Pocessing Formation %d of the current ECLIPSE data:\n' , G);
                Geometry = mcomputeGeometry(grdecl_processed(G));
                temp = obj.ComputeGeometry(Geometry);
                CornerPointGridData(n).Nx              = temp.Nx;
                CornerPointGridData(n).Ny              = temp.Ny;
                CornerPointGridData(n).Nz              = temp.Nz;
                CornerPointGridData(n).N_TotalCells    = temp.N_TotalCells;
                CornerPointGridData(n).N_ActiveCells   = temp.N_ActiveCells;
                CornerPointGridData(n).N_Nodes         = temp.N_Nodes;
                CornerPointGridData(n).N_InternalFaces = temp.N_InternalFaces;
                CornerPointGridData(n).N_ExternalFaces = temp.N_ExternalFaces;
                CornerPointGridData(n).Nodes           = temp.Nodes;
                CornerPointGridData(n).Cells           = temp.Cells;
                CornerPointGridData(n).Internal_Faces  = temp.Internal_Faces;
                CornerPointGridData(n).External_Faces  = temp.External_Faces;

                %% SECTION 2.1: Porosity Data
                if any(strcmp(fieldnames(grdecl_raw),'PORO')) && ~isempty(grdecl_raw.PORO)
                    CornerPointGridData(n).Porosity = grdecl_raw.PORO(Geometry.cells.indexMap);
                end
                
                %% SECTION 2.2: Permeability Data
                if any(strcmp(fieldnames(grdecl_raw),'PERMX')) && ~isempty(grdecl_raw.PERMX)
                    if ~strcmp(fieldnames(grdecl_raw),'PERMY')
                        grdecl_raw.PERMY = grdecl_raw.PERMX;
                    end
                    if ~strcmp(fieldnames(grdecl_raw),'PERMZ')
                        grdecl_raw.PERMZ = grdecl_raw.PERMY;
                    end
                    perm = [grdecl_raw.PERMX, grdecl_raw.PERMY, grdecl_raw.PERMZ];
                    perm = perm(Geometry.cells.indexMap,:);
                    % Convert K values from to Darcy to Metric (if not already in metric)
                    CornerPointGridData(n).Permeability = perm .* milli * darcy;
                    CornerPointGridData(n).PermUnit = 'm2';
                    CornerPointGridData(n).PermScale = 'Linear';
                end
            end
        end
        function CornerPointGridData = ComputeGeometry(obj,Geometry)
            %% SECTION 1.1: CELLS NODES (X, Y, Z COORDINATES)
            fprintf('---> Pocessing Nodes ... ');
            CornerPointGridData.Nodes = Geometry.nodes.coords;
            fprintf('Done!\n');
            
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
            FtC_2 = obj.ObtainFaceIndices(faces, pos, cells);                                                              % c = cells
            
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
            CtC_2 = obj.ObtainCellNeighbors(faces, pos, cells);                                                         % c = cells
            CtC_2 = sort(CtC_2,2);
            
            fprintf('---> Pocessing Cells ... ');
            Cell.Vertices    = NtC;
            Cell.Centroid    = Geometry.cells.centroids;
            Cell.Volume      = Geometry.cells.volumes;
            Cell.Faces       = num2cell(FtC_2,2);
            Cell.Faces       = cellfun( @(x) unique(x,'stable') , Cell.Faces , 'UniformOutput' , false);
            Cell.N_Faces     = cellfun( @length , Cell.Faces );
            Cell.Neighbors   = num2cell(CtC_2,2);
            Cell.Neighbors   = cellfun( @(x) unique(x,'stable') , Cell.Neighbors , 'UniformOutput' , false);
            Cell.N_Neighbors =  cellfun( @length , Cell.Neighbors );
            
            Cell.dx = zeros(size(Cell.Volume));
            Cell.dy = zeros(size(Cell.Volume));
            Cell.dz = zeros(size(Cell.Volume));
            for i = 1 : size(Cell.Volume,1)
                Cell.dx(i) = max( CornerPointGridData.Nodes( Cell.Vertices(i,:),1) ) - min( CornerPointGridData.Nodes( Cell.Vertices(i,:),1) );
                Cell.dy(i) = max( CornerPointGridData.Nodes( Cell.Vertices(i,:),2) ) - min( CornerPointGridData.Nodes( Cell.Vertices(i,:),2) );
                Cell.dz(i) = max( CornerPointGridData.Nodes( Cell.Vertices(i,:),3) ) - min( CornerPointGridData.Nodes( Cell.Vertices(i,:),3) );
            end
            
            CornerPointGridData.Cells = Cell;
            fprintf('Done!\n');
            
            %% SECTION 1.3: INTERNAL FACES (FACES CONNECTED TO CELLS)
            nodes = Geometry.faces.nodes;
            pos   = Geometry.faces.nodePos;
            faces = 1:Geometry.faces.num ;
            NtF_3 = obj.ObtainNodeIndices(nodes, pos, faces);                                                             % f = Faces
            
            NF = linspace(1, Geometry.faces.num, Geometry.faces.num)';                                            % Create Face Index Vector (Total Number of Faces)
            
            % Assembly Matrix with Face Geometry Data
            IF = [NF, Geometry.faces.neighbors, Geometry.faces.areas, Geometry.faces.centroids, Geometry.faces.normals NtF_3];
            LI = (IF(:,2) == 0)|(IF(:,3) == 0);                                                                   % Delete External Faces
            IF(LI,:) = [];                                                                                        % Delete External Faces
            % Create Centroid Vector: Face Centroid - Cell Centroid
            c_vec = [(Geometry.faces.centroids(IF(:,1),:) - Geometry.cells.centroids(IF(:,2),:)) (Geometry.faces.centroids(IF(:,1),:) - Geometry.cells.centroids(IF(:,3),:))];
            
            % Internal Faces Data: Face Index + Face Area + Face Centroid + Face Normal + Cell Neighbor + Centroid Vector
            fprintf('---> Pocessing Internal Faces ... ');
            Internal_Face.FullIndex          = IF(:,1);
            Internal_Face.Area               = IF(:,4);
            Internal_Face.Centroid           = IF(:,5:7);
            Internal_Face.Nvec               = IF(:,8:10);
            Internal_Face.CellNeighbor1Index = IF(:,2);
            Internal_Face.CellNeighbor1Vec   = c_vec(:,1:3);
            Internal_Face.CellNeighbor2Index = IF(:,3);
            Internal_Face.CellNeighbor2Vec   = c_vec(:,4:6);
            Internal_Face.Vertices           = num2cell( IF(:,11:end) , 2 );
            Internal_Face.Vertices           = cellfun( @(x) unique(x,'stable') , Internal_Face.Vertices , 'UniformOutput' , false);
            Internal_Face.N_Vertices         = cellfun( @length , Internal_Face.Vertices );
            CornerPointGridData.Internal_Faces = Internal_Face;
            fprintf('Done!\n');
            
            %% SECTION 1.4: EXTERNAL FACES (FACES AT THE EXTERNAL BOUNDARIES OF THE GRID)
            % Assembly Matrix with Face Geometry Data
            EF = [NF, Geometry.faces.areas, Geometry.faces.centroids, Geometry.faces.normals, Geometry.faces.neighbors NtF_3];
            LI = (EF(:,9) ~= 0)&(EF(:,10) ~= 0);                                                                 % Delete Internal Faces
            EF(LI,:) = [];                                                                                       % Delete Internal Faces
            NC = EF(:,9)+ EF(:,10);                                                                              % Delete Cell Neighboors  == 0
            
            % External Faces Data: Face Index + Face Area + Face Centroid + Face Normal + Cell Neighbor + Centroid Vector
            fprintf('---> Pocessing External Faces ... ');
            External_Face.FullIndex         = EF(:,1);
            External_Face.Area              = EF(:,2);
            External_Face.Centroid          = EF(:,3:5);
            External_Face.Nvec              = EF(:,6:8);
            External_Face.CellNeighborIndex = NC;                   % An external face has only one connection (to only one cell).
            External_Face.CellNeighborVec   = Geometry.faces.centroids(EF(:,1),:) - Geometry.cells.centroids(NC);
            External_Face.Vertices          = num2cell( EF(:,11:end) , 2 );
            External_Face.Vertices          = cellfun( @(x) unique(x,'stable') , External_Face.Vertices , 'UniformOutput' , false);
            External_Face.N_Vertices        = cellfun( @length , External_Face.Vertices );
            CornerPointGridData.External_Faces = External_Face;
            fprintf('Done!\n');
            
            % Adding extra statistical data
            CornerPointGridData.Nx = Geometry.cartDims(1);
            CornerPointGridData.Ny = Geometry.cartDims(2);
            CornerPointGridData.Nz = Geometry.cartDims(3);
            CornerPointGridData.N_TotalCells = prod(Geometry.cartDims);
            CornerPointGridData.N_ActiveCells = Geometry.cells.num;
            CornerPointGridData.N_InactiveCells = CornerPointGridData.N_TotalCells - CornerPointGridData.N_ActiveCells;
            CornerPointGridData.N_TotalFaces = Geometry.faces.num;
            CornerPointGridData.N_InternalFaces = length(CornerPointGridData.Internal_Faces.Area);
            CornerPointGridData.N_ExternalFaces = length(CornerPointGridData.External_Faces.Area);
            CornerPointGridData.N_Nodes = Geometry.nodes.num;
            fprintf('\n');
        end
        function CornerPointGridData = ReadfromTXT(obj,CornerPointGridRockPropertiesFile)
            fid = fopen(obj.File, 'r');
            GeometryMatrix = textscan(fid, '%s', 'Delimiter', '\n');
            GeometryMatrix = GeometryMatrix{1};
            fclose(fid);
            if nargin > 1 && ~isempty(CornerPointGridRockPropertiesFile)
                fid = fopen(CornerPointGridRockPropertiesFile, 'r');
                RockPropertiesMatrix = textscan(fid, '%s', 'Delimiter', '\n');
                RockPropertiesMatrix = RockPropertiesMatrix{1};
                fclose(fid);
            else
                RockPropertiesMatrix = [];
            end
            
            fprintf('Reading the CornerPointGrid input file:\n');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'RESERVOIR_GRID_NX');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "RESERVOIR_GRID_NX" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.Nx = str2double( GeometryMatrix{index+1} );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'RESERVOIR_GRID_NY');
            if isempty(index)
                error('The keyword "RESERVOIR_GRID_NY" is missing. Please check the CornerPointGrid input file!\n');
            end
            index = find(~cellfun('isempty', temp));
            CornerPointGridData.Ny = str2double( GeometryMatrix{index+1} );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'RESERVOIR_GRID_NZ');
            if isempty(index)
                error('The keyword "RESERVOIR_GRID_NZ" is missing. Please check the CornerPointGrid input file!\n');
            end
            index = find(~cellfun('isempty', temp));
            CornerPointGridData.Nz = str2double( GeometryMatrix{index+1} );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            index = find(strcmp(GeometryMatrix, 'ACTIVE_CELLS'));
            if isempty(index)
                error('The keyword "ACTIVE_CELLS" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.N_ActiveCells = str2double( GeometryMatrix{index+1} );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            index = find(strcmp(GeometryMatrix, 'INACTIVE_CELLS'));
            if isempty(index)
                error('The keyword "ACTIVE_CELLS" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.N_InctiveCells = str2double( GeometryMatrix{index+1} );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            index = find(strcmp(GeometryMatrix, 'TOTAL_CELLS'));
            if isempty(index)
                error('The keyword "ACTIVE_CELLS" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.N_TotalCells = str2double( GeometryMatrix{index+1} );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'N_INTERNAL_FACES');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "N_INTERNAL_FACES" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.N_InternalFaces = str2double( GeometryMatrix{index+1} );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'N_EXTERNAL_FACES');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "N_EXTERNAL_FACES" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.N_ExternalFaces = str2double( GeometryMatrix{index+1} );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'N_NODES');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "N_NODES" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.N_Nodes = str2double( GeometryMatrix{index+1} );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'MAX_NUMBER_OF_FACES_PER_CELL');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "MAX_NUMBER_OF_FACES_PER_CELL" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.MaxNumFacesPerCell = str2double( GeometryMatrix{index+1} );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'MAX_NUMBER_OF_NEIGHBORS_PER_CELL');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "MAX_NUMBER_OF_NEIGHBORS_PER_CELL" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.MaxNumNeighborsPerCell = str2double( GeometryMatrix{index+1} );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'MAX_NUMBER_OF_NODES_PER_FACE');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "MAX_NUMBER_OF_NODES_PER_FACE" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.MaxNumVertxPerFace = str2double( GeometryMatrix{index+1} );
            
            %%% Reading nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'NODES_COORDINATES');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "NODES_COORDINATES" is missing. Please check the CornerPointGrid input file!\n');
            end
            
            fprintf('---> Reading Nodes ... ');
            TEMP = GeometryMatrix(index+2:index+2+CornerPointGridData.N_Nodes-1);
            splitStr = regexp(TEMP, ',', 'split');
            splitStr = str2double( vertcat( splitStr{:} ) );
            Nodes = splitStr(:,2:4);
            fprintf('Done!\n');
            CornerPointGridData.Nodes = Nodes;
            
            %%% Reading cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'CELL_GEOMETRY');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "CELL_GEOMETRY" is missing. Please check the CornerPointGrid input file!\n');
            end

            mnf = CornerPointGridData.MaxNumFacesPerCell;
            mnn = CornerPointGridData.MaxNumNeighborsPerCell;
            fprintf('---> Reading Cells ... ');
            TEMP = GeometryMatrix(index+2:index+2+CornerPointGridData.N_ActiveCells-1);
            splitStr = regexp(TEMP, ',', 'split');
            splitStr = str2double( vertcat( splitStr{:} ) );
            Cell.Vertices = splitStr(:,2:9);
            Cell.Centroid = splitStr(:,10:12);
            Cell.Volume = splitStr(:,13);
            Cell.Faces = num2cell( splitStr(:,14:14+mnf-1) , 2 );
            Cell.Faces = cellfun( @(x) unique(x,'stable') , Cell.Faces , 'UniformOutput' , false);
            Cell.N_Faces = cellfun( @length , Cell.Faces );
            Cell.Neighbors = num2cell( splitStr(:,14+mnf:14+mnf+mnn-1) , 2 );
            Cell.Neighbors = cellfun( @(x) unique(x,'stable') , Cell.Neighbors , 'UniformOutput' , false);
            Cell.N_Neighbors =  cellfun( @length , Cell.Neighbors );
            
            Cell.dx = zeros(size(Cell.Volume));
            Cell.dy = zeros(size(Cell.Volume));
            Cell.dz = zeros(size(Cell.Volume));
            for i = 1 : obj.CornerPointGridData.N_ActiveCells
                Cell.dx(i) = max( CornerPointGridData.Nodes( Cell.Vertices(i,:),1) ) - min( CornerPointGridData.Nodes( Cell.Vertices(i,:),1) );
                Cell.dy(i) = max( CornerPointGridData.Nodes( Cell.Vertices(i,:),2) ) - min( CornerPointGridData.Nodes( Cell.Vertices(i,:),2) );
                Cell.dz(i) = max( CornerPointGridData.Nodes( Cell.Vertices(i,:),3) ) - min( CornerPointGridData.Nodes( Cell.Vertices(i,:),3) );
            end
            
            fprintf('Done!\n');
            CornerPointGridData.Cells = Cell;
            
            %%% Reading internal faces  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'INTERNAL_FACE_GEOMETRY');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "INTERNAL_FACE_GEOMETRY" is missing. Please check the CornerPointGrid input file!\n');
            end
            
            mnv = CornerPointGridData.MaxNumVertxPerFace;
            fprintf('---> Reading Internal Faces ... ');
            TEMP = GeometryMatrix(index+2:index+2+CornerPointGridData.N_InternalFaces-1);
            splitStr = regexp(TEMP, ',', 'split');
            splitStr = str2double( vertcat( splitStr{:} ) );
            Internal_Face.FullIndex          = splitStr(:,1);
            Internal_Face.Area               = splitStr(:,2);
            Internal_Face.Centroid           = splitStr(:,3:5);
            Internal_Face.Nvec               = splitStr(:,6:8);
            Internal_Face.CellNeighbor1Index = splitStr(:,9);
            Internal_Face.CellNeighbor1Vec   = splitStr(:,10:12);
            Internal_Face.CellNeighbor2Index = splitStr(:,13);
            Internal_Face.CellNeighbor2Vec   = splitStr(:,14:16);
            Internal_Face.Vertices = num2cell( splitStr(:,17:17+mnv-1) , 2);
            Internal_Face.Vertices = cellfun( @(x) unique(x,'stable') , Internal_Face.Vertices , 'UniformOutput' , false);
            Internal_Face.N_Vertices = cellfun( @length , Internal_Face.Vertices );
            fprintf('Done!\n');
            CornerPointGridData.Internal_Faces = Internal_Face;
            
            %%% Reading external faces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(GeometryMatrix, 'EXTERNAL_FACE_GEOMETRY');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "EXTERNAL_FACE_GEOMETRY" is missing. Please check the CornerPointGrid input file!\n');
            end
            
            mnv = CornerPointGridData.MaxNumVertxPerFace;
            fprintf('---> Reading External Faces ... ');
            TEMP = GeometryMatrix(index+2:index+2+CornerPointGridData.N_ExternalFaces-1);
            splitStr = regexp(TEMP, ',', 'split');
            splitStr = str2double( vertcat( splitStr{:} ) );
            External_Face.FullIndex         = splitStr(:,1);
            External_Face.Area              = splitStr(:,2);
            External_Face.Centroid          = splitStr(:,3:5);
            External_Face.Nvec              = splitStr(:,6:8);
            External_Face.CellNeighborIndex = splitStr(:,9);     % An external face has only one connection (to only one cell).
            External_Face.CellNeighborVec   = splitStr(:,10:12);
            External_Face.Vertices = num2cell( splitStr(:,13:13+mnv-1) , 2);
            External_Face.Vertices = cellfun( @(x) unique(x,'stable') , External_Face.Vertices , 'UniformOutput' , false);
            External_Face.N_Vertices = cellfun( @length , External_Face.Vertices );
            fprintf('Done!\n');
            CornerPointGridData.External_Faces = External_Face;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Reading the rock properties (porosity and permeability)
            
            if isempty(RockPropertiesMatrix)
                return;
            end
            fprintf('Reading the CornerPointGrid rock properties data input file:\n');
            
            temp = strfind(RockPropertiesMatrix, 'RESERVOIR_GRID_NX');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "RESERVOIR_GRID_NX" is missing. Please check the CornerPointGrid rock properties input file!\n');
            end
            if CornerPointGridData.Nx ~= str2double( RockPropertiesMatrix{index+1} )
            	error('The number of grids in x (Nx) of the CornerPointGrid rock properties input file does not match with CornerPointGrid data file!\n');
            end
            
            temp = strfind(RockPropertiesMatrix, 'RESERVOIR_GRID_NY');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "RESERVOIR_GRID_NY" is missing. Please check the CornerPointGrid input file!\n');
            end
            if CornerPointGridData.Ny ~= str2double( RockPropertiesMatrix{index+1} )
            	error('The number of grids in y (Ny) of the CornerPointGrid rock properties input file does not match with CornerPointGrid data file!\n');
            end
                
            temp = strfind(RockPropertiesMatrix, 'RESERVOIR_GRID_NZ');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "RESERVOIR_GRID_NZ" is missing. Please check the CornerPointGrid input file!\n');
            end
            if CornerPointGridData.Nz ~= str2double( RockPropertiesMatrix{index+1} )
            	error('The number of grids in z (Nz) of the CornerPointGrid rock properties input file does not match with CornerPointGrid data file!\n');
            end
                
            temp = strfind(RockPropertiesMatrix, 'ACTIVE_CELLS');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "ACTIVE_CELLS" is missing. Please check the CornerPointGrid input file!\n');
            end
            if CornerPointGridData.N_ActiveCells ~= str2double( RockPropertiesMatrix{index+1} )
            	error('The number of ative cells of the CornerPointGrid rock properties input file does not match with CornerPointGrid data file!\n');
            end
            
            %%% Reading porosity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(RockPropertiesMatrix, 'POROSITY_DATA');
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                fprintf('---> Reading Porosity ... ');
                TEMP = RockPropertiesMatrix(index+2:index+2+CornerPointGridData.N_ActiveCells-1);
                splitStr = regexp(TEMP,',','split');
                splitStr = vertcat(splitStr{:});
                CornerPointGridData.Porosity = str2double(splitStr(:,end));
                fprintf('Done!\n');
            end
            
            %%% Reading permeability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            temp = strfind(RockPropertiesMatrix, 'PERMEABILITY_DATA');
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                fprintf('---> Reading Permeability ... ');
                TEMP = RockPropertiesMatrix(index+2:index+2+CornerPointGridData.N_ActiveCells-1);
                splitStr = regexp(TEMP,',','split');
                splitStr = vertcat(splitStr{:});
                CornerPointGridData.Permeability = str2double(splitStr(:,2:end));
                fprintf('Done!\n');
                
                temp = strfind(RockPropertiesMatrix, 'PERMEABILITY_UNIT');
                index = find(~cellfun('isempty', temp));
                if ~isempty(index)
                    CornerPointGridData.PermUnit = RockPropertiesMatrix{index+1};
                else
                    warning('The keyword "PERMEABILITY_UNIT" is missing. SI unit "m2" is set by default.\n');
                    CornerPointGridData.PermUnit = 'm2';
                end
                 
                temp = strfind(RockPropertiesMatrix, 'PERMEABILITY_SCALE');
                index = find(~cellfun('isempty', temp));
                if ~isempty(index)
                    CornerPointGridData.PermScale = RockPropertiesMatrix{index+1};
                else
                    warning('The keyword "PERMEABILITY_SCALE" is missing. Linear scale is set by default.\n');
                    CornerPointGridData.PermScale = 'Linear';
                end
                 
                if ~strcmp(CornerPointGridData.PermScale, 'Linear') && ~strcmp(CornerPointGridData.PermScale, 'Logarithmic')
                    error('The permeability scale is either Linear or Logarithmic. Please check the input file.\n');
                end
            end
        end
        function [f, present] = ObtainNodeIndices(obj, nodes, pos, faces)
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
            tmp         = isfinite(f);
            nnode       = sum(tmp,1);
            ind         = sub2ind(size(f),nnode,1:size(f,2));
            tmp         = repmat(f(ind),size(f,1),1);
            f(isnan(f)) = tmp(isnan(f));
            
            f = f .';
        end
        function [c, present] = ObtainFaceIndices(obj, faces, pos, cells)
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
            tmp         = isfinite(c);
            nnode       = sum(tmp,1);
            ind         = sub2ind(size(c),nnode,1:size(c,2));
            tmp         = repmat(c(ind),size(c,1),1);
            c(isnan(c)) = tmp(isnan(c));
            
            c = c .';
        end
        function [c, present] = ObtainCellNeighbors(obj, faces, pos, cells)
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
            tmp         = isfinite(c);
            nnode       = sum(tmp,1);
            ind         = sub2ind(size(c),nnode,1:size(c,2));
            tmp         = repmat(c(ind),size(c,1),1);
            c(isnan(c)) = tmp(isnan(c));
            
            c = c .';
        end
    end
end