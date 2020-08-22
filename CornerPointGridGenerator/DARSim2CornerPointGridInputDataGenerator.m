%% CornerPointGrid Converter for DARSim2
% Based on MRST (Matlab Reservoir Simulation Toolbox), it reads the input file(in an Eclipse file format)and generates the input file 
% with the cell geometry that can be read by DARSim2.
% Author: Janio Piguave
% Modified on: 2020/08/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The functions "processGRDECL" and "computeGeometry" can be substituted by "mprocessGRDECL" and "mcomputeGeometry"
% respectively, in order to use C-accelerated library which speeds up the process of ECLISPE models data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DARSim2CornerPointGridInputDataGenerator(Directory,FileName)
%% CORNER POINT GRID GEOMETRY DATA: INPUT FILE GENERATION
% Read subset of ECLIPSE GRID file. Return: CartDims + ZCORN + COORD + ACTNUM
Reader = reader_eclipse(Directory, FileName);
CornerPointGridData = Reader.ReadGRDECL('all');

% Some models have more than one disconnected formations. Use "length(Geometries)" as
% the maximum counter for the following for-loop to process all of them.
for n = 1:length(CornerPointGridData)
    G = CornerPointGridData(n).Index;
    
    %% SECTION 1: OUTPUT TXT FILE FOR GRID GEOMETRY DATA
    disp( '******************* Writing the data into output text files *********************' );
    OutputFileName = 'CornerPointGrid_DARSim_InputData';
    disp(['Writing into file ', OutputFileName, ' #', num2str(G)]);
    
    % Writing the header
    fid = fopen(strcat(Directory,'\',OutputFileName,'_',FileName,'_',num2str(G),'.txt') , 'w+' );
    fprintf(fid, '%% *******************************************************\n');
    fprintf(fid, '%% CORNER POINT GRID DATA FOR DARSIM 2 RESERVOIR SIMULATOR\n');
    fprintf(fid, '%% *******************************************************\n');
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'RESERVOIR_GRID_NX\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).Nx );
    fprintf(fid, '\n');
    fprintf(fid, 'RESERVOIR_GRID_NY\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).Ny );
    fprintf(fid, '\n');
    fprintf(fid, 'RESERVOIR_GRID_NZ\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).Nz );
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'TOTAL_CELLS\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).N_TotalCells );
    fprintf(fid, '\n');
    fprintf(fid, 'ACTIVE_CELLS\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).N_ActiveCells );
    fprintf(fid, '\n');
    fprintf(fid, 'INACTIVE_CELLS\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).N_InactiveCells );
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'N_FACES\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).N_TotalFaces );
    fprintf(fid, '\n');
    fprintf(fid, 'N_INTERNAL_FACES\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).N_InternalFaces );
    fprintf(fid, '\n');
    fprintf(fid, 'N_EXTERNAL_FACES\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).N_ExternalFaces );
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'N_NODES\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).N_Nodes );
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'MAX_NUMBER_OF_FACES_PER_CELL\n');
    fprintf(fid, '%d\n', max(CornerPointGridData(n).Cells.N_Faces) );
    fprintf(fid, '\n');
    fprintf(fid, 'MAX_NUMBER_OF_NEIGHBORS_PER_CELL\n');
    fprintf(fid, '%d\n', max(CornerPointGridData(n).Cells.N_Neighbors) );
    fprintf(fid, '\n');
    fprintf(fid, 'MAX_NUMBER_OF_NODES_PER_FACE\n');
    fprintf(fid, '%d\n', max(CornerPointGridData(n).Internal_Faces.N_Vertices) );
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, '\n');
    
    % Writing the nodes
    TXT_Nodes = horzcat( [1:CornerPointGridData(n).N_Nodes]' , CornerPointGridData(n).Nodes);
    fprintf(fid, '%% Section 1: Nodes Coordinates\n');
    fprintf(fid, '%% [Node Index] + [Nodes Coordinates (x;y;z)]\n');
    fprintf(fid, '\n');
    fprintf(fid, 'NODES_COORDINATES\n');
    fprintf(fid,'%s %7s %16s %13s\n','Node No.  ','x','y','z');
    fprintf(fid,'%7d,    %6f , %6f , %6f\n', TXT_Nodes' );
    fprintf(fid, '\n');
    fprintf(fid, '----------------------------------------------------------------------------------------------------------\n');
    
    % Writing the cells
    maxNumFaces = max(CornerPointGridData(n).Cells.N_Faces);
    CellsFaces = cell2mat( cellfun(@(x) [x ones(1,maxNumFaces-numel(x))*x(end)] , CornerPointGridData(n).Cells.Faces , 'uni' , 0 ) );
    
    maxNumNeighbors = max(CornerPointGridData(n).Cells.N_Neighbors);
    CellsNeighbors = cell2mat( cellfun(@(x) [x ones(1,maxNumNeighbors-numel(x))*x(end)] , CornerPointGridData(n).Cells.Neighbors , 'uni' , 0 ) );
    
    TXT_Cells = horzcat( [1:CornerPointGridData(n).N_ActiveCells]' , CornerPointGridData(n).Cells.Vertices , CornerPointGridData(n).Cells.Centroid , ...
                         CornerPointGridData(n).Cells.Volume , CellsFaces , CellsNeighbors );
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
    for f = 1 : maxNumFaces
        FormatSpec = strcat(FormatSpec,"%6d,");
    end
    FormatSpec = char(FormatSpec);
    FormatSpec(end)=[];
    FormatSpec = strcat(FormatSpec,"    ,    ");
    for f = 1 : maxNumNeighbors
        FormatSpec = strcat(FormatSpec,"%6d,");
    end
    FormatSpec = char(FormatSpec);
    FormatSpec(end)=[];
    FormatSpec = strcat(FormatSpec,'\n');
    
    fprintf(fid,FormatSpec,TXT_Cells');
    fprintf(fid, '\n');
    fprintf(fid, '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
    
    % Writing the internal faces
    maxNumVertices = max(CornerPointGridData(n).Internal_Faces.N_Vertices);
    FacesVertices = cell2mat( cellfun(@(x) [x ones(1,maxNumVertices-numel(x))*x(end)] , CornerPointGridData(n).Internal_Faces.Vertices , 'uni' , 0 ) );
    TXT_Internal_Faces = horzcat( CornerPointGridData(n).Internal_Faces.FullIndex , CornerPointGridData(n).Internal_Faces.Area , ...
                         CornerPointGridData(n).Internal_Faces.Centroid           , CornerPointGridData(n).Internal_Faces.Nvec , ...
                         CornerPointGridData(n).Internal_Faces.CellNeighbor1Index , CornerPointGridData(n).Internal_Faces.CellNeighbor1Vec , ...
                         CornerPointGridData(n).Internal_Faces.CellNeighbor2Index , CornerPointGridData(n).Internal_Faces.CellNeighbor2Vec , ...
                         FacesVertices);
    fprintf(fid, '%% Section 3: Internal Faces (shared with two cells)\n');
    fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell 1] + [Centroid Vector 1(x,y,z)] + [Neighboring Cell 2] + [Centroid Vector 2(x,y,z)] + [Node Indices]\n');
    fprintf(fid, '%% NC: Neighbor Cell\n');
    fprintf(fid, '\n');
    fprintf(fid, 'INTERNAL_FACE_GEOMETRY\n');
    fprintf(fid,'%s %18s %39s %47s %22s %33s %18s %33s %27s\n',...
                'Face No.','Face Area','Face Centroid(x;y;z)','Face Normal(x;y;z)','NC1','Centroid Vector1(x:y:z)','NC2','Centroid Vector2(x,y,z)','Node Indices');
    
    FormatSpec = "%7.0d    ,    %12.6f    ,    %12.6f,%12.6f,%12.6f    ,   %12.6f,%12.6f,%12.6f    ,   %6.0d,%12.6f,%12.6f,%12.6f    ,   %6.0d,%12.6f,%12.6f,%12.6f    ,   ";
    for f = 1 : maxNumVertices
        FormatSpec = strcat(FormatSpec,"%6d,");
    end
    FormatSpec = char(FormatSpec);
    FormatSpec(end)=[];
    FormatSpec = strcat(FormatSpec,'\n');
    
    fprintf(fid,FormatSpec,TXT_Internal_Faces');
    fprintf(fid, '\n');
    fprintf(fid, '------------------------------------------------------------------------------------------------------------------------------------------------\n');
    
    % Writing the external faces
    maxNumVertices = max(CornerPointGridData(n).External_Faces.N_Vertices);
    FacesVertices = cell2mat( cellfun(@(x) [x ones(1,maxNumVertices-numel(x))*x(end)] , CornerPointGridData(n).External_Faces.Vertices , 'uni' , 0 ) );
    TXT_External_Faces = horzcat( CornerPointGridData(n).External_Faces.FullIndex , CornerPointGridData(n).External_Faces.Area , ...
                         CornerPointGridData(n).External_Faces.Centroid           , CornerPointGridData(n).External_Faces.Nvec , ...
                         CornerPointGridData(n).External_Faces.CellNeighborIndex  , CornerPointGridData(n).External_Faces.CellNeighborVec , ...
                         FacesVertices);
    fprintf(fid, '%% Section 4: External Faces (at boundaries | no shared with cells)\n');
    fprintf(fid, '%% [Face Index] + [Face Area] + [Face Centroid (x,y,z)] + [Face Normal(x,y,z)] + [Neighboring Cell] + [Centroid Vector (x,y,z)] + [Face Topology]\n');
    fprintf(fid, '%% NC: Neighbor Cell\n');
    fprintf(fid, '\n');
    fprintf(fid, 'EXTERNAL_FACE_GEOMETRY\n');
    fprintf(fid,'%s %18s %39s %47s %22s %33s %29s\n',...
                'Face No.','Face Area','Face Centroid(x;y;z)','Face Normal(x;y;z)','NC','Centroid Vector(x,y,z)','Face Topology');
    
    FormatSpec = "%7.0d    ,    %12.6f    ,    %12.6f,%12.6f,%12.6f    ,   %12.6f,%12.6f,%12.6f    ,   %6.0d,%12.6f,%12.6f,%12.6f    ,   ";
    for f = 1 : maxNumVertices
        FormatSpec = strcat(FormatSpec,"%6d,");
    end
    FormatSpec = char(FormatSpec);
    FormatSpec(end)=[];
    FormatSpec = strcat(FormatSpec,'\n');
    
    fprintf(fid,FormatSpec,TXT_External_Faces');
    
    fclose(fid);
    
    %% SECTION 2: OUTPUT VTK FILE FOR PARAVIEW
    OutputFileName = 'CornerPointGrid_DARSim_VTK_Plot';
    disp(['Writing into file ', OutputFileName, ' #', num2str(G)]);
    fid = fopen(strcat(Directory,'\',OutputFileName,'_',FileName,'_',num2str(G),'.vtk') , 'w+' );
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'DARSim2 Reservoir Simulator\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    
    % Nodes: Set Depth from shallow to deep
    VTK_Nodes = [CornerPointGridData(n).Nodes(:,1:2) , max(CornerPointGridData(n).Nodes(:,3))-CornerPointGridData(n).Nodes(:,3)];
    fprintf(fid, ['POINTS ' num2str(CornerPointGridData(n).N_Nodes) ' double\n']);
    fprintf(fid, '%f %f %f\n', VTK_Nodes');
    fprintf(fid, '\n');
    
    % Cells: Indexing is from 0 to N-1 , 8 nodes number per cell + vertices 
    VTK_Cells = horzcat(8*ones(CornerPointGridData(n).N_ActiveCells,1) , CornerPointGridData(n).Cells.Vertices-1 );
    fprintf(fid, ['CELLS ' num2str(CornerPointGridData(n).N_ActiveCells) ' ' num2str(CornerPointGridData(n).N_ActiveCells*(1+8)) '\n']);
    fprintf(fid, '%d %d %d %d %d %d %d %d %d\n',VTK_Cells');
    fprintf(fid, '\n');
    
    fprintf(fid, ['CELL_TYPES ' num2str(CornerPointGridData(n).N_ActiveCells) '\n']);
    fprintf(fid, '%d ', 11 * ones(1,CornerPointGridData(n).N_ActiveCells) );
    fprintf(fid, '\n\n');
    
    % Writing porosity & permeability into VTK file
    fprintf(fid, 'CELL_DATA %d\n', CornerPointGridData(n).N_ActiveCells);
    fprintf(fid, ' \n');
    fprintf(fid, strcat('SCALARS Porosity double 1\n'));
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid,'%1.5e ', CornerPointGridData(n).Porosity');
    fprintf(fid, '\n\n');
    fprintf(fid, strcat('SCALARS PermX double 1\n'));
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid,'%1.5e ', CornerPointGridData(n).Permeability(:,1)');
    fprintf(fid, '\n\n');
    fprintf(fid, strcat('SCALARS PermY double 1\n'));
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid,'%1.5e ', CornerPointGridData(n).Permeability(:,2)');
    fprintf(fid, '\n\n');
    fprintf(fid, strcat('SCALARS PermZ double 1\n'));
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid,'%1.5e ', CornerPointGridData(n).Permeability(:,3)');
    fprintf(fid, '\n\n');
    fclose(fid);
    
    %% SECTION 3: OUTPUT TXT FILE FOR ROCK PROPERTIES
    OutputFileName = 'CornerPointGrid_DARSim_RockPropertiesData';
    disp(['Writing into file ', OutputFileName, ' #', num2str(G)]);
    fid = fopen(strcat(Directory,'\',OutputFileName,'_',FileName,'_',num2str(G),'.txt') , 'w+' );
    fprintf(fid, '%% Rock Properties of the Norne Field\n');
    fprintf(fid, '%% Rock Properties Values of Active Cells of the Corner Point Grid Model\n');
    fprintf(fid, '%% Permeability tensor is assumed to be diagonal\n');
    fprintf(fid, '%% Vertical permeability (Kz) equaling one-tenth of the horizontal permeability (Kx, Ky)\n');
    fprintf(fid, '%% Only the x-component Kx is given in the data file. Kx = Ky. Kz = 0.1 Kx\n');
    fprintf(fid, '\n');
    fprintf(fid, '%% The Grid Resolution of the Reservoir:\n');
    fprintf(fid, 'RESERVOIR_GRID_NX\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).Nx );
    fprintf(fid, 'RESERVOIR_GRID_NY\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).Ny );
    fprintf(fid, 'RESERVOIR_GRID_NZ\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).Nz );
    fprintf(fid, 'ACTIVE_CELLS\n');
    fprintf(fid, '%d\n', CornerPointGridData(n).N_ActiveCells );
    
    % Writing the porosity data
    if any(strcmp(fieldnames(CornerPointGridData(n)),'Porosity')) && ~isempty(CornerPointGridData(n).Porosity)
        TXT_Porosity = horzcat( [1:CornerPointGridData(n).N_ActiveCells]' , CornerPointGridData(n).Porosity);
        fprintf(fid, '\n');
        fprintf(fid, '%% Section 1: Porosity\n');
        fprintf(fid, '%% [Cell Index] + [Porosity]\n');
        fprintf(fid, '%% Units: Fraction\n');
        fprintf(fid, 'POROSITY_DATA\n');
        fprintf(fid,'%s %s\n', 'Cell No.  ','Porosity');
        fprintf(fid,'%6.0d   ,  %.4f\n', TXT_Porosity' );
    end
    
    % Writing the permeability data
    if any(strcmp(fieldnames(CornerPointGridData(n)),'Permeability')) && ~isempty(CornerPointGridData(n).Permeability)
        TXT_Permeability = horzcat( [1:CornerPointGridData(n).N_ActiveCells]' , CornerPointGridData(n).Permeability);
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
        fprintf(fid,'%6.0d   ,  %.4e,%.4e,%.4e\n', TXT_Permeability');
    end
    fclose(fid);    
end
end