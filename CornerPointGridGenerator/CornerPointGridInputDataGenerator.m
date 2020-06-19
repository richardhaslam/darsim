%% CornerPointGrid Converter for DARSim2
% Based on MRST (Matlab Reservoir Simulation Toolbox), it reads the input file(in an Eclipse file format)
% and generates the input file with the cell geometry that can be read by DARSim2.
% 
% Author: Janio Piguave
% Modified on: 2020/01/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CornerPointGridInputDataGenerator(Directory,FileName)
%% CORNER POINT GRID GEOMETRY DATA: INPUT FILE GENERATION
% Read subset of ECLIPSE GRID file. Return: CartDims + ZCORN + COORD + ACTNUM
InputFile = strcat(Directory,'\',FileName);
grdecl = readGRDECL(InputFile);
[x, y, z] = buildCornerPtNodes(grdecl);                                    % Extract Corner Point Nodes of cells (8 nodes * X,Y,Z(coordinates))
Geometries = processGRDECL(grdecl, 'Verbose', true);                                % Compute grid topology and geometry from pillar grid description
% Compute geometry of  cells: centroids, volumes, areas
% grdecl.ACTNUM(Geometries(2).cells.indexMap)= 0;                                     % Norne Field: Delete Small Section
% G = computeGeometry(G(1));                                                 % Norne Field: Compute geometry of  cells: centroids, volumes, areas

for G = 1 : length(Geometries)
    Geometry = computeGeometry(Geometries(G));
    ACTNUM = grdecl.ACTNUM;
    % Removing the active cells of other geometries (if any) from the raw data (grdecl.ACTNUM)
    for i = [1:G-1,G+1:length(Geometries)]
        ACTNUM(Geometries(i).cells.indexMap)= 0;
    end
    
    % Reshape the cell data based on the format requirements of DARSim2 input file:
    X = inputdataDARSim(x);
    Y = inputdataDARSim(y);
    Z = inputdataDARSim(z);
    
    % SECTION 1: CELLS NODES (X, Y, Z) + CENTROIDS + VECTORS A, B, C (X, Y, Z COORDINATES)
    % Create Cell Index Vector (Number of Cells)
    NC = linspace(1, Geometry.cartDims(1)*Geometry.cartDims(2)*Geometry.cartDims(3), Geometry.cartDims(1)*Geometry.cartDims(2)*Geometry.cartDims(3))';
    
    Cell_Nodes = [NC, double(ACTNUM),...
                  X(:,2), Y(:,2), Z(:,2), X(:,6), Y(:,6), Z(:,6),...
                  X(:,1), Y(:,1), Z(:,1), X(:,5), Y(:,5), Z(:,5),...
                  X(:,4), Y(:,4), Z(:,4), X(:,8), Y(:,8), Z(:,8),...
                  X(:,3), Y(:,3), Z(:,3), X(:,7), Y(:,7), Z(:,7)];
    
    % Only Active Cells (based on ACTNUM data)
    LI = Cell_Nodes(:,2) == 0;                            % Logical Index
    Cell_Nodes(LI,:) = [];                                % Delete Cells that are not active
    Cell_Nodes(:,2) = [];                                 % Delete Columns of Active Cells
    Cell_Data = [Cell_Nodes, Geometry.cells.centroids, Geometry.cells.volumes];
    
    % SECTION 2: INTERNAL FACES (FACES CONNECTED TO CELLS)
    NF = linspace(1, Geometry.faces.num, Geometry.faces.num)';                               % Create Face Index Vector (Total Number of Faces)
    % Assembly Matrix with Face Geometry Data
    IF = [NF, Geometry.faces.neighbors, Geometry.faces.areas, Geometry.faces.centroids, Geometry.faces.normals];
    LI = (IF(:,2) == 0)|(IF(:,3) == 0);                                        % Delete External Faces
    IF(LI,:) = [];                                                             % Delete External Faces
    % Create Centroid Vector: Face Centroid - Cell Centroid
    c_vec = [Geometry.faces.centroids(IF(:,1),:) - Geometry.cells.centroids(IF(:,2),:), Geometry.faces.centroids(IF(:,1),:) - Geometry.cells.centroids(IF(:,3),:)];
    % Internal Faces Data: Face Index + Face Area + Face Centroid + Face Normal + Cell Neighbor + Centroid Vector
    IF2 = [IF(:,1) IF(:,4:10) IF(:,2) c_vec(:,1:3) IF(:,3) c_vec(:,4:6)];
    
    % SECTION 3: EXTERNAL FACES (FACES AT THE EXTERNAL BOUNDARIES OF THE GRID)
    % Assembly Matrix with Face Geometry Data
    EF = [NF, Geometry.faces.areas, Geometry.faces.centroids, Geometry.faces.normals, Geometry.faces.neighbors];
    LI = (EF(:,9) ~= 0)&(EF(:,10) ~= 0);                                       % Delete Internal Faces
    EF(LI,:) = [];                                                             % Delete Internal Faces
    EF2 = [EF(:,1:8) (EF(:,9)+EF(:,10))];                                      % Delete Cell Neighboors  == 0
    % External Faces Data: Face Index + Face Area + Face Centroid + Face Normal + Cell Neighbor + Centroid Vector
    EF3 = [EF2 Geometry.faces.centroids(EF2(:,1),:) - Geometry.cells.centroids(EF2(:,9))];
    %% CORNER POINT GRID ROCK PROPERTIES DATA: INPUT FILE GENERATION
    % POROSITY
    if sum(strcmp(fieldnames(grdecl),'PORO'))
        poro = grdecl.PORO(Geometry.cells.indexMap);                                % SAIGUP Model | Norne Field
        poro_text = [Geometry.cells.indexMap poro];
    else
        poro = load('Porosity.txt')';                                      % Johansen Formation
        poro = poro(Geometry.cells.indexMap);                                     % Johansen Formation
    end
    
    % PERMEABILITY
    if sum(strcmp(fieldnames(grdecl),'PERMX')) % such as in "SAIGU"P and "NOrne"
        if sum(strcmp(fieldnames(grdecl),'PERMZ'))
            perm = [grdecl.PERMX, grdecl.PERMY, grdecl.PERMZ];
        else
            perm = bsxfun(@times, [1 1 1], grdecl.PERMX);
        end
    else % such as in "Johansen"
        perm = load('Permeability.txt')';
        perm = bsxfun(@times, [1 1 1], perm);
    end
    perm = perm(Geometry.cells.indexMap,:);
    
    perm = perm .* milli * darcy;                                              % Convert K values from to Darcy to Metric
    perm_txt = [Geometry.cells.indexMap perm];
    
    rock = makeRock(Geometry, perm, poro);                                            % Create K Tensor
    [K, i, j] = permTensor(rock, Geometry.griddim);                                   % Create K Tensor
    K_text = [Geometry.cells.indexMap K];                                             % Create K Tensor
    
    %% OUTOUT FILE 1: GRID GEOMETRY DATA
    OutputFileName = 'CornerPointGrid_DARSim_InputData';
    disp( '******************* Writing the data into output text file *********************' );
    disp(['Writing into file ', OutputFileName, '#', G]);
    
    fid = fopen(strcat(Directory,OutputFileName,'_',num2str(G),'.txt') , 'w+' );
    fprintf(fid, '%% Node Locations for each one of the cells, Nx * Ny *Nz\n');
    fprintf(fid, '%% Coordinates X,Y,Z of the eight nones for each one of the cells\n');
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
    %% OUTPUT FILE 2: ROCK PROPERTIES DATA
    OutputFileName = 'CornerPointGrid_DARSim_RockPropertiesData';
    disp( '******************* Writing the data into output text file *********************' );
    disp(['Writing into file ', OutputFileName, '#', G]);
    
    fid = fopen(strcat(Directory, OutputFileName, '_', num2str(G),'.txt'), 'w+');
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

end

%%--------------------------------------------------------------------------
function A = inputdataDARSim(B)
C = size(B);
D = reshape(permute(reshape(B, C(1), 2, [], C(3)), [1,3,2,4]), [], 2, C(3));
A = []; 
E=1;
for k=2:2:size(D,3)
    F(:,:,E)=[D(:,:,k-1),D(:,:,k)];                             
    E=E+1;
end
for k=1:size(D,3)/2
    A=[A;F(:,:,k)];                                               
end
A=A';
A=reshape(A,8,[])';
end
%--------------------------------------------------------------------------
%% End of main function
end