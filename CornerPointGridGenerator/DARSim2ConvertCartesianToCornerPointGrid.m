%% ConvertCartesianToCornerPointGrid for DARSim
% This script receives the geometry input from a Cartesian domain and
% creates CornerPointGrid input data file (txt) for DARSim simulator.
% 
% Author: Mousa HosseiniMehr
% Modified on: 2020/04/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DARSim2ConvertCartesianToCornerPointGrid(Directory, InputFileName)
%% Reading the DARSim Main Input File
InputFile = strcat(Directory,'/',InputFileName);
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
cartGeometry = cartGrid([Nx, Ny, Nz], [Lx, Ly, Lz]);                                         % Create Cartesian Grid
Geometry = mcomputeGeometry(cartGeometry);
Reader = reader_eclipse(InputFile);
CornerPointGridData = Reader.ComputeGeometry(Geometry);                                  % Compute Geomtery of the Cartesian Grid

%% SECTION 1: OUTPUT TXT FILE FOR GRID GEOMETRY DATA
    disp( '******************* Writing the data into output text files *********************' );
    OutputFileName = 'CornerPointGrid_DARSim_InputData';
    disp(['Writing into file ', OutputFileName, ' #', num2str(G)]);
    
    % Writing the header
    fid = fopen(strcat(Directory,'\',OutputFileName,'_',FileName,'.txt') , 'w+' );
    fprintf(fid, '%% *******************************************************\n');
    fprintf(fid, '%% CORNER POINT GRID DATA FOR DARSIM 2 RESERVOIR SIMULATOR\n');
    fprintf(fid, '%% *******************************************************\n');
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'RESERVOIR_GRID_NX\n');
    fprintf(fid, '%d\n', CornerPointGridData.Nx );
    fprintf(fid, '\n');
    fprintf(fid, 'RESERVOIR_GRID_NY\n');
    fprintf(fid, '%d\n', CornerPointGridData.Ny );
    fprintf(fid, '\n');
    fprintf(fid, 'RESERVOIR_GRID_NZ\n');
    fprintf(fid, '%d\n', CornerPointGridData.Nz );
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'TOTAL_CELLS\n');
    fprintf(fid, '%d\n', CornerPointGridData.N_TotalCells );
    fprintf(fid, '\n');
    fprintf(fid, 'ACTIVE_CELLS\n');
    fprintf(fid, '%d\n', CornerPointGridData.N_ActiveCells );
    fprintf(fid, '\n');
    fprintf(fid, 'INACTIVE_CELLS\n');
    fprintf(fid, '%d\n', CornerPointGridData.N_InactiveCells );
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'N_FACES\n');
    fprintf(fid, '%d\n', CornerPointGridData.N_TotalFaces );
    fprintf(fid, '\n');
    fprintf(fid, 'N_INTERNAL_FACES\n');
    fprintf(fid, '%d\n', CornerPointGridData.N_InternalFaces );
    fprintf(fid, '\n');
    fprintf(fid, 'N_EXTERNAL_FACES\n');
    fprintf(fid, '%d\n', CornerPointGridData.N_ExternalFaces );
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'N_NODES\n');
    fprintf(fid, '%d\n', CornerPointGridData.N_Nodes );
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'MAX_NUMBER_OF_FACES_PER_CELL\n');
    fprintf(fid, '%d\n', max(CornerPointGridData.Cells.N_Faces) );
    fprintf(fid, '\n');
    fprintf(fid, 'MAX_NUMBER_OF_NEIGHBORS_PER_CELL\n');
    fprintf(fid, '%d\n', max(CornerPointGridData.Cells.N_Neighbors) );
    fprintf(fid, '\n');
    fprintf(fid, 'MAX_NUMBER_OF_NODES_PER_FACE\n');
    fprintf(fid, '%d\n', max(CornerPointGridData.Internal_Faces.N_Vertices) );
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, '\n');
    
    % Writing the nodes
    TXT_Nodes = horzcat( [1:CornerPointGridData.N_Nodes]' , CornerPointGridData.Nodes);
    fprintf(fid, '%% Section 1: Nodes Coordinates\n');
    fprintf(fid, '%% [Node Index] + [Nodes Coordinates (x;y;z)]\n');
    fprintf(fid, '\n');
    fprintf(fid, 'NODES_COORDINATES\n');
    fprintf(fid,'%s %7s %16s %13s\n','Node No.  ','x','y','z');
    fprintf(fid,'%7d,    %6f , %6f , %6f\n', TXT_Nodes' );
    fprintf(fid, '\n');
    fprintf(fid, '----------------------------------------------------------------------------------------------------------\n');
    
    % Writing the cells
    maxNumFaces = max(CornerPointGridData.Cells.N_Faces);
    CellsFaces = cell2mat( cellfun(@(x) [x ones(1,maxNumFaces-numel(x))*x(end)] , CornerPointGridData.Cells.Faces , 'uni' , 0 ) );
    
    maxNumNeighbors = max(CornerPointGridData.Cells.N_Neighbors);
    CellsNeighbors = cell2mat( cellfun(@(x) [x ones(1,maxNumNeighbors-numel(x))*x(end)] , CornerPointGridData.Cells.Neighbors , 'uni' , 0 ) );
    
    TXT_Cells = horzcat( [1:CornerPointGridData.N_ActiveCells]' , CornerPointGridData.Cells.Vertices , CornerPointGridData.Cells.Centroid , ...
                         CornerPointGridData.Cells.Volume , CellsFaces , CellsNeighbors );
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
    maxNumVertices = max(CornerPointGridData.Internal_Faces.N_Vertices);
    FacesVertices = cell2mat( cellfun(@(x) [x ones(1,maxNumVertices-numel(x))*x(end)] , CornerPointGridData.Internal_Faces.Vertices , 'uni' , 0 ) );
    TXT_Internal_Faces = horzcat( CornerPointGridData.Internal_Faces.FullIndex , CornerPointGridData.Internal_Faces.Area , ...
                         CornerPointGridData.Internal_Faces.Centroid           , CornerPointGridData.Internal_Faces.Nvec , ...
                         CornerPointGridData.Internal_Faces.CellNeighbor1Index , CornerPointGridData.Internal_Faces.CellNeighbor1Vec , ...
                         CornerPointGridData.Internal_Faces.CellNeighbor2Index , CornerPointGridData.Internal_Faces.CellNeighbor2Vec , ...
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
    maxNumVertices = max(CornerPointGridData.External_Faces.N_Vertices);
    FacesVertices = cell2mat( cellfun(@(x) [x ones(1,maxNumVertices-numel(x))*x(end)] , CornerPointGridData.External_Faces.Vertices , 'uni' , 0 ) );
    TXT_External_Faces = horzcat( CornerPointGridData.External_Faces.FullIndex , CornerPointGridData.External_Faces.Area , ...
                         CornerPointGridData.External_Faces.Centroid           , CornerPointGridData.External_Faces.Nvec , ...
                         CornerPointGridData.External_Faces.CellNeighborIndex  , CornerPointGridData.External_Faces.CellNeighborVec , ...
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
    fid = fopen(strcat(Directory,'\',OutputFileName,'_',FileName,'.vtk') , 'w+' );
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'DARSim2 Reservoir Simulator\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    
    % Nodes: Set Depth from shallow to deep
    VTK_Nodes = [CornerPointGridData.Nodes(:,1:2) , max(CornerPointGridData.Nodes(:,3))-CornerPointGridData.Nodes(:,3)];
    fprintf(fid, ['POINTS ' num2str(CornerPointGridData.N_Nodes) ' double\n']);
    fprintf(fid, '%f %f %f\n', VTK_Nodes');
    fprintf(fid, '\n');
    
    % Cells: Indexing is from 0 to N-1 , 8 nodes number per cell + vertices 
    VTK_Cells = horzcat(8*ones(CornerPointGridData.N_ActiveCells,1) , CornerPointGridData.Cells.Vertices-1 );
    fprintf(fid, ['CELLS ' num2str(CornerPointGridData.N_ActiveCells) ' ' num2str(CornerPointGridData.N_ActiveCells*(1+8)) '\n']);
    fprintf(fid, '%d %d %d %d %d %d %d %d %d\n',VTK_Cells');
    fprintf(fid, '\n');
    
    fprintf(fid, ['CELL_TYPES ' num2str(CornerPointGridData.N_ActiveCells) '\n']);
    fprintf(fid, '%d ', 11 * ones(1,CornerPointGridData.N_ActiveCells) );
    fprintf(fid, '\n\n');
    
    % Writing porosity & permeability into VTK file
    fprintf(fid, 'CELL_DATA %d\n', CornerPointGridData.N_ActiveCells);
    fprintf(fid, ' \n');
    fprintf(fid, strcat('SCALARS Porosity double 1\n'));
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid,'%1.5e ', CornerPointGridData.Porosity');
    fprintf(fid, '\n\n');
    fprintf(fid, strcat('SCALARS PermX double 1\n'));
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid,'%1.5e ', CornerPointGridData.Permeability(:,1)');
    fprintf(fid, '\n\n');
    fprintf(fid, strcat('SCALARS PermY double 1\n'));
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid,'%1.5e ', CornerPointGridData.Permeability(:,2)');
    fprintf(fid, '\n\n');
    fprintf(fid, strcat('SCALARS PermZ double 1\n'));
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid,'%1.5e ', CornerPointGridData.Permeability(:,3)');
    fprintf(fid, '\n\n');
    fclose(fid);
end