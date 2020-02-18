% Reader for Matteo's input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 December 2017
%Last modified: 14 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reader_darsim2 < reader
    properties
        InputMatrix
        SettingsMatrix
        FractureMatrix
        CornerPointGridMatrix
        CornerPointGridRockPropertiesMatrix
    end
    methods
        function obj = reader_darsim2(dir, file, permdirectory)
            obj@reader(dir, file, permdirectory);          
        end
        function ReadInputFile(obj, Builder)
            %% Read Input File
            fileID = fopen(obj.File, 'r');
            % Read lines from input file
            matrix = textscan(fileID, '%s', 'Delimiter', '\n');
            obj.InputMatrix = matrix{1};
            fclose(fileID);
            % Remove lines which are commented (contain --)
            Commented = startsWith(obj.InputMatrix, '--');
            obj.InputMatrix(Commented) = {'--'}; % removing the string if it is commented.
            
            %% Read Settings File
            SettingsFile = strcat(obj.Directory, '/SimulatorSettings.txt');
            fileID = fopen(SettingsFile, 'r');
            matrix = textscan(fileID, '%s', 'Delimiter', '\n');
            obj.SettingsMatrix = matrix{1};
            fclose(fileID);
            % Remove lines which are commented (contain --)
            Commented = startsWith(obj.SettingsMatrix, '--');
            obj.SettingsMatrix(Commented) = {'--'}; % removing the string if it is commented.
            
            %% Read Fracture File
            fractured = sum(contains(obj.InputMatrix, 'FRACTURED'));
            if fractured 
                FractureFile = strcat(obj.Directory, '/Fracture_Output.txt');
                fileID = fopen(FractureFile, 'r');
				fprintf('Reading fracture file ...');
                matrix = textscan(fileID, '%s', 'Delimiter', '\n');
                obj.FractureMatrix = matrix{1};
                fclose(fileID);
				fprintf('--> Successful.\n\n');
            end
            [Builder.SimulationInput, Builder.SimulatorSettings] = ReadInformation(obj);
        end
        function [SimulationInput, SimulatorSettings] = ReadInformation(obj)
            %% INPUT FILE
            % 1. Name of the Problem 
            temp = strfind(obj.InputMatrix, 'TITLE'); % Search a specific string and find all rows containing matches
            SimulationInput.ProblemName = char(obj.InputMatrix(find(~cellfun('isempty', temp)) + 1));
			
            % 2. Total time
            temp = strfind(obj.InputMatrix, 'TOTALTIME');
            xv = find(~cellfun('isempty', temp));
            if isempty(xv)
                error('The keyword "TOTALTIME" is missing. Please check the input file!\n');
            end
            Time = strsplit( obj.InputMatrix{xv + 1} , ' ');
            Hours = 0; Minutes = 0; Seconds = 0;
            Days = str2double(Time{1});
            if length(Time)>1
                if strcmp(Time{2},'DAYS')
                    if strcmp(Time{4},'HOURS')  ,  Hours   = str2double(Time{3});  end
                    if strcmp(Time{6},'MINUTES'),  Minutes = str2double(Time{5});  end
                    if strcmp(Time{8},'SECONDS'),  Seconds = str2double(Time{7});  end
                end
            end
            SimulationInput.TotalTime = Days*24*3600 + Hours*3600 + Minutes*60 + Seconds;
			
            % 3. Read reservoir properties
            SimulationInput.ReservoirProperties = obj.ReadReservoirProperties();
            
            % 4. Read fluid properties
            SimulationInput.FluidProperties = obj.ReadFluidProperties();
                                  
            % 5. Initial conditions
            temp = strfind(obj.InputMatrix, 'INIT');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "INIT" is missing. Please check the input file!\n');
            end
            SimulationInput.Init = str2double(strsplit(char(obj.InputMatrix(index + 1))));
            
            % 6. Read wells info
            SimulationInput.WellsInfo = obj.ReadWellsInfo(SimulationInput);
            
            % 7. Read Fractures' Properties (if any)
            temp = strfind(obj.InputMatrix, 'FRACTURED'); % Check the main input file and look for FRACTURED
            Fractured = find(~cellfun('isempty', temp));
            if isempty(Fractured)
                SimulationInput.FracturesProperties.Fractured = 0;
                SimulationInput.FracturesProperties.NrOfFrac = 0;
            else
                 SimulationInput.FracturesProperties = obj.ReadFracturesProperties();
            end
            
            %% SIMULATOR'S SETTINGS
            SimulatorSettings = obj.ReadSimulatorSettings(SimulationInput);
        end
        function ReservoirProperties = ReadReservoirProperties(obj)
            %%%%%%%%%%%%%PROPERTIES OF THE RESERVOIR%%%%%%%%%%%%%%%%
            temp = strfind(obj.InputMatrix, 'GRID_DISCRETIZATION');
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                ReservoirProperties.Discretization = obj.InputMatrix{index+1};
            else
                warning('The keyword "GRID_DISCRETIZATION" is missing. Cartesian discretization is set by default.\n');
                ReservoirProperties.Discretization = 'Cartesian';
            end
            
            if strcmp(ReservoirProperties.Discretization,'CornerPointGrid')
                % Read the cornerpointgrid file
                FileName = obj.InputMatrix{index+2};
                if isempty(FileName)
                    FileName = 'CornerPointGrid_DARSim_InputData.txt';
                end
                CornerPointGridFile = strcat(obj.Directory, '/', FileName);
                fileID = fopen(CornerPointGridFile, 'r');
                matrix = textscan(fileID, '%s', 'Delimiter', '\n');
                obj.CornerPointGridMatrix = matrix{1};
                fclose(fileID);
                
                if isfile(strcat(obj.Directory, '/','CornerPointGridData.mat'))
                    fprintf('"CornerPointGridData.mat" file already exists. No need to load the CornerPointGrid data input file.\n');
                    load(strcat(obj.Directory, '/','CornerPointGridData.mat'),'ReservoirProperties');
                else
                    ReservoirProperties.CornerPointGridData = obj.ReadCornerPointGridData();
                    save(strcat(obj.Directory, '/','CornerPointGridData.mat'),'ReservoirProperties');
                end
                
                ReservoirProperties.Grid.N(1) = ReservoirProperties.CornerPointGridData.Nx;
                ReservoirProperties.Grid.N(2) = ReservoirProperties.CornerPointGridData.Ny;
                ReservoirProperties.Grid.N(3) = ReservoirProperties.CornerPointGridData.Nz;
                ReservoirProperties.Grid.N_ActiveCells = ReservoirProperties.CornerPointGridData.N_ActiveCells;
                
                % For now temporarily, we get Lx,Ly,LZ of the reservoir
                % from the main input file. Soon, we will read this from
                % the CornerPointGrid data.
                temp = strfind(obj.InputMatrix, 'DIMENS'); % Search a specific string and find all rows containing matches
                index = find(~cellfun('isempty', temp));
                if isempty(index)
                    error('The keyword "DIMENS" is missing. Please check the input file!\n');
                end
                ReservoirProperties.size = [str2double(obj.InputMatrix{index+1});...
                                            str2double(obj.InputMatrix{index+2});
                                            str2double(obj.InputMatrix{index+3})];
                
            elseif strcmp(ReservoirProperties.Discretization,'Cartesian')
                % Assume it is Cartesian
                % 1. size of the reservoir
                temp = strfind(obj.InputMatrix, 'DIMENS'); % Search a specific string and find all rows containing matches
                index = find(~cellfun('isempty', temp));
                if isempty(index)
                    error('The keyword "DIMENS" is missing. Please check the input file!\n');
                end
                ReservoirProperties.size = [str2double(obj.InputMatrix{index+1});...
                                            str2double(obj.InputMatrix{index+2});
                                            str2double(obj.InputMatrix{index+3})];
                % 2. GridSize
                temp = strfind(obj.InputMatrix, 'SPECGRID');
                index = find(~cellfun('isempty', temp));
                if isempty(index)
                    error('The keyword "SPECGRID" is missing. Please check the input file!\n');
                end
                ReservoirProperties.Grid.N = [str2double(obj.InputMatrix{index+1});...
                                              str2double(obj.InputMatrix{index+2});
                                              str2double(obj.InputMatrix{index+3})];
            else
                error('The discretization method should either be "Cartesian" or "CornerPointGrid". Check the input file!\n');
            end
            
            
            % 3. Permeability
            if strcmp(ReservoirProperties.Discretization,'CornerPointGrid')
                temp = strfind(obj.InputMatrix, 'ROCKPROPERTIES_FILE');
                index = find(~cellfun('isempty', temp));
                FileName = obj.InputMatrix{index+1};
                if isempty(FileName)
                    FileName = 'CornerPointGrid_DARSim_RockPropertiesData.txt';
                end
                CornerPointGridRockPropertiesFile = strcat(obj.Directory, '/', FileName);
                fileID = fopen(CornerPointGridRockPropertiesFile, 'r');
                matrix = textscan(fileID, '%s', 'Delimiter', '\n');
                obj.CornerPointGridRockPropertiesMatrix = matrix{1};
                fclose(fileID);
                
                if isfile(strcat(obj.Directory, '/','CornerPointGridRockPropertiesData.mat'))
                    fprintf('"CornerPointGridRockPropertiesData.mat" file already exists. No need to load the CornerPointGridRockProperties data input file.\n');
                    load(strcat(obj.Directory, '/','CornerPointGridRockPropertiesData.mat'),'ReservoirProperties');
                else
                    ReservoirProperties.CornerPointGridRockPropertiesData = obj.ReadCornerPointGridRockPropertiesData();
                    save(strcat(obj.Directory, '/','CornerPointGridRockPropertiesData.mat'),'ReservoirProperties');
                end
                
            else
                
                strcmp(ReservoirProperties.Discretization,'Cartesian')
                temp = strfind(obj.InputMatrix, 'PERMEABILITY_UNIT');
                index = find(~cellfun('isempty', temp));
                if ~isempty(index)
                    ReservoirProperties.PermUnit = obj.InputMatrix{index+1};
                else
                    warning('The keyword "PERMEABILITY_UNIT" is missing. SI unit "m2" is set by default.\n');
                    ReservoirProperties.PermUnit = 'm2';
                end
                 
                temp = strfind(obj.InputMatrix, 'PERMEABILITY_SCALE');
                index = find(~cellfun('isempty', temp));
                if ~isempty(index)
                    ReservoirProperties.PermScale = obj.InputMatrix{index+1};
                else
                    warning('The keyword "PERMEABILITY_SCALE" is missing. Linear scale is set by default.\n');
                    ReservoirProperties.PermScale = 'Linear';
                end
                 
                if ~strcmp(ReservoirProperties.PermScale, 'Linear') && ~strcmp(ReservoirProperties.PermScale, 'Logarithmic')
                    error('The permeability scale is either Linear or Logarithmic. Please check the input file.\n');
                end
                 
                temp = strfind(obj.InputMatrix, 'PERMEABILITY_CONTRAST_REDUCTION');
                index = find(~cellfun('isempty', temp));
                if ~isempty(index)
                    ReservoirProperties.PermContrastReduction = 1;
                    temp1 = strfind(obj.InputMatrix, 'PERMEABILITY_CONTRAST_MEAN');
                    index1 = find(~cellfun('isempty', temp1));
                    if ~isempty(index1)
                        ReservoirProperties.PermContrastMean = obj.InputMatrix{index1+1};
                    else
                        ReservoirProperties.PermContrastMean = 'Default';
                    end
                    temp1 = strfind(obj.InputMatrix, 'PERMEABILITY_CONTRAST_ORDER');
                    index1 = find(~cellfun('isempty', temp1));
                    if ~isempty(index1)
                        ReservoirProperties.PermContrastOrder = str2double(obj.InputMatrix{index1+1});
                    else
                        ReservoirProperties.PermContrastOrder = 3;
                    end
                else
                    ReservoirProperties.PermContrastReduction = 0;
                end
			
                perm = zeros(3, 1);
                temp = strfind(obj.InputMatrix, 'PERMX');
                perm(1) = find(~cellfun('isempty', temp));
                if isempty(perm(1))
                    error('The keyword "PERMX" is missing. Please check the input file!\n');
                end
                temp = strfind(obj.InputMatrix, 'PERMY');
                perm(2) = find(~cellfun('isempty', temp));
                if isempty(perm(2))
                    error('The keyword "PERMY" is missing. Please check the input file!\n');
                end
                temp = strfind(obj.InputMatrix, 'PERMZ');
                perm(3) = find(~cellfun('isempty', temp));
                if isempty(perm(3))
                    error('The keyword "PERMX" is missing. Please check the input file!\n');
                end
                DimChar = {'x' , 'y' , 'z'};
                for i=1:3
                    if contains(obj.InputMatrix(perm(i) - 1), 'INCLUDE')
                        ReservoirProperties.PermInclude(i) = 1;
                        ReservoirProperties.PermFile{i} = strcat(obj.PermDirectory, char(obj.InputMatrix(perm(i) +1)));
                    else
                        ReservoirProperties.PermInclude(i) = 0;
                        ReservoirProperties.Perm(i) = str2double(obj.InputMatrix(perm(i) +1));
                        if isnan(ReservoirProperties.Perm(i))
                            error('The permeability value in %s direction is not valid. Please check the input file.\n',DimChar{i});
                        end
                    end
                end
			
                %%%% Homogenized/upscaled permeabilities for each coarsening level %%%%%%%
                temp = strfind(obj.SettingsMatrix, 'DLGR');
                x = find(~cellfun('isempty', temp));
                if isempty(temp)
                    DLGR = 0;
                else
                    DLGR = str2double(obj.SettingsMatrix(x+1));
                end
                if DLGR
                    temp = strfind(obj.SettingsMatrix, 'ADM');
                    adm = find(~cellfun('isempty', temp));
                    if str2double(obj.SettingsMatrix(adm + 1)) == 1
                        temp = strfind(obj.SettingsMatrix, 'LEVELS');
                        x = find(~cellfun('isempty', temp));
                        CoarseningLevel = str2double(obj.SettingsMatrix(x+1));
                        for L = 1:CoarseningLevel
                            temp = strfind( obj.InputMatrix , strcat('PERM_COARSE_L',num2str(L),'_X') );
                            permx = find(~cellfun('isempty', temp));
                            temp = strfind( obj.InputMatrix , strcat('PERM_COARSE_L',num2str(L),'_Y') );
                            permy = find(~cellfun('isempty', temp));
                            temp = strfind( obj.InputMatrix , strcat('PERM_COARSE_L',num2str(L),'_Z') );
                            permz = find(~cellfun('isempty', temp));
                            if ~isempty(permx)
                                ReservoirProperties.CoarsePermFile{L,1} = strcat(obj.PermDirectory, char(obj.InputMatrix(permx+1)));
                            end
                            if ~isempty(permy)
                                ReservoirProperties.CoarsePermFile{L,2} = strcat(obj.PermDirectory, char(obj.InputMatrix(permy+1)));
                            end
                            if ~isempty(permz)
                                ReservoirProperties.CoarsePermFile{L,3} = strcat(obj.PermDirectory, char(obj.InputMatrix(permz+1)));
                            end
                        end
                    end
                end
                
                % 4. Porosity
                temp = strfind(obj.InputMatrix, 'POR');
                index = find(~cellfun('isempty', temp));
                ReservoirProperties.phi = str2double(obj.InputMatrix(index + 1));
                
            end
            
            
            % 5. Temperature
            temp = strfind(obj.InputMatrix, 'TEMPERATURE (K)');
            index = find(~cellfun('isempty', temp));
            ReservoirProperties.Temperature = str2double(obj.InputMatrix(index(1) + 1));
            
            % Reading the number of phases
            temp = strfind(obj.InputMatrix, 'FLUID MODEL');
            index = find(~cellfun('isempty', temp));
            NofPhases = str2double(obj.InputMatrix{index+3});
            
            % 6. Rock Compressibility
            temp = strfind(obj.InputMatrix, 'COMPRESSIBILITY');
            index_comp = find(~cellfun('isempty', temp));
            ReservoirProperties.Compressibility = str2double(obj.InputMatrix{index_comp+((NofPhases+1)*2)});
             if isnan(ReservoirProperties.Compressibility)
                 ReservoirProperties.Compressibility = 0; % Default value if not defined [1/Pa]
             end
            
            % 7. Rock Density
            temp = strfind(obj.InputMatrix, 'DENSITY');
            index_density = find(~cellfun('isempty', temp));
            ReservoirProperties.RockDensity = str2double(obj.InputMatrix{index_density+((NofPhases+1)*2)});
            if isnan(ReservoirProperties.RockDensity)
                ReservoirProperties.RockDensity = 2750; % Default value if not defined [kg/m3]
            end
            
            % 8. Rock Conductivity
            temp = strfind(obj.InputMatrix, 'CONDUCTIVITY');
            index_conduc = find(~cellfun('isempty', temp));
            if isempty(index_conduc)
                ReservoirProperties.RockConductivity = 4; % Default value if not defined [W/m/K]
            else
                ReservoirProperties.RockConductivity = str2double(obj.InputMatrix{index_conduc+((NofPhases+1)*2)});
            end
            
            % 9. Rock Specific Heat
            temp = strfind(obj.InputMatrix, 'SPECIFIC HEAT');
            index_spec_heat = find(~cellfun('isempty', temp));
            if isempty(index_spec_heat)
                ReservoirProperties.SpecificHeat = 790; % Default value if not defined [J/Kg/K]
            else
                ReservoirProperties.SpecificHeat = str2double(obj.InputMatrix{index_spec_heat+(Permeability(NofPhases+1)*2)});
            end
            
        end
        function FracturesProperties = ReadFracturesProperties(obj)
            %%%%%%%%%%%%%PROPERTIES OF THE FRACTURES%%%%%%%%%%%%%%%%
            FracturesProperties.Fractured = 1;
            temp = strfind(obj.FractureMatrix, 'NUM_FRACS');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "NUM_FRACS" is missing. Please check the Fracture_Output file!\n');
            end
            temp = strsplit(obj.FractureMatrix{index},' ');
            FracturesProperties.NrOfFrac = str2double( temp{end} );
            
        end
        function CornerPointGridData = ReadCornerPointGridData(obj)
            fprintf('Reading the CornerPointGrid input file:\n');
            temp = strfind(obj.CornerPointGridMatrix, 'RESERVOIR_GRID_NX');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "RESERVOIR_GRID_NX" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.Nx = str2double( obj.CornerPointGridMatrix{index+1} );
            
            temp = strfind(obj.CornerPointGridMatrix, 'RESERVOIR_GRID_NY');
            if isempty(index)
                error('The keyword "RESERVOIR_GRID_NY" is missing. Please check the CornerPointGrid input file!\n');
            end
            index = find(~cellfun('isempty', temp));
            CornerPointGridData.Ny = str2double( obj.CornerPointGridMatrix{index+1} );
            
            temp = strfind(obj.CornerPointGridMatrix, 'RESERVOIR_GRID_NZ');
            if isempty(index)
                error('The keyword "RESERVOIR_GRID_NZ" is missing. Please check the CornerPointGrid input file!\n');
            end
            index = find(~cellfun('isempty', temp));
            CornerPointGridData.Nz = str2double( obj.CornerPointGridMatrix{index+1} );
            
            temp = strfind(obj.CornerPointGridMatrix, 'ACTIVE_CELLS');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "ACTIVE_CELLS" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.N_ActiveCells = str2double( obj.CornerPointGridMatrix{index+1} );
            
            temp = strfind(obj.CornerPointGridMatrix, 'N_INTERNAL_FACES');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "N_INTERNAL_FACES" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.N_InternalFaces = str2double( obj.CornerPointGridMatrix{index+1} );
            
            temp = strfind(obj.CornerPointGridMatrix, 'N_EXTERNAL_FACES');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "N_EXTERNAL_FACES" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridData.N_ExternalFaces = str2double( obj.CornerPointGridMatrix{index+1} );
            
            % Reading internal face information line by line
            temp = strfind(obj.CornerPointGridMatrix, 'INTERNAL_FACE_GEOMETRY');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "INTERNAL_FACE_GEOMETRY" is missing. Please check the CornerPointGrid input file!\n');
            end
            
            Internal_Face.Area = zeros(CornerPointGridData.N_InternalFaces , 1);
            Internal_Face.Centroid = zeros(CornerPointGridData.N_InternalFaces , 3);
            Internal_Face.Nvec = zeros(CornerPointGridData.N_InternalFaces , 3);
            Internal_Face.CellNeighbor1Index = zeros(CornerPointGridData.N_InternalFaces , 1);
            Internal_Face.CellNeighbor1Vec = zeros(CornerPointGridData.N_InternalFaces , 3);
            Internal_Face.CellNeighbor2Index = zeros(CornerPointGridData.N_InternalFaces , 1);
            Internal_Face.CellNeighbor2Vec = zeros(CornerPointGridData.N_InternalFaces , 3);
            
            fprintf('---> Internal Face ');
            for i = 1 : CornerPointGridData.N_InternalFaces
                if (i>1),  fprintf(repmat('\b', 1, 13));  end
                fprintf('%06d/%06d',i,CornerPointGridData.N_InternalFaces);
                Temp = strsplit(obj.CornerPointGridMatrix{index+2-1+i},{' , '});
                Internal_Face.Area(i) =  str2double( Temp{2} );
                Internal_Face.Centroid(i,:) = str2double( strsplit( Temp{3} , ';' ) );
                Internal_Face.Nvec(i,:) = str2double( strsplit( Temp{4} , ';' ) );
                Internal_Face.CellNeighbor1Index(i) = str2double( Temp{5} );
                Internal_Face.CellNeighbor1Vec(i,:) = str2double( strsplit( Temp{6} , ';' ) );
                Internal_Face.CellNeighbor2Index(i) = str2double( Temp{7} );
                Internal_Face.CellNeighbor2Vec(i,:) = str2double( strsplit( Temp{8} , ';' ) );
            end
            fprintf('\n');
            CornerPointGridData.Internal_Face = Internal_Face;
            
            % Reading external face information line by line
            temp = strfind(obj.CornerPointGridMatrix, 'EXTERNAL_FACE_GEOMETRY');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "EXTERNAL_FACE_GEOMETRY" is missing. Please check the CornerPointGrid input file!\n');
            end
            External_Face.Area = zeros(CornerPointGridData.N_ExternalFaces , 1);
            External_Face.Centroid = zeros(CornerPointGridData.N_ExternalFaces , 3);
            External_Face.Nvec = zeros(CornerPointGridData.N_ExternalFaces , 3);
            External_Face.CellNeighborIndex = zeros(CornerPointGridData.N_ExternalFaces , 1);
            External_Face.CellNeighborVec = zeros(CornerPointGridData.N_ExternalFaces , 3);
            
            fprintf('---> External Face ');
            for i = 1 : CornerPointGridData.N_ExternalFaces
                if (i>1),  fprintf(repmat('\b', 1, 13));  end
                fprintf('%06d/%06d',i,CornerPointGridData.N_ExternalFaces);
                Temp = strsplit(obj.CornerPointGridMatrix{index+2-1+i},{' , '});
                External_Face.Area(i) =  str2double( Temp{2} );
                External_Face.Centroid(i,:) = str2double( strsplit( Temp{3} , ';' ) );
                External_Face.Nvec(i,:) = str2double( strsplit( Temp{4} , ';' ) );
                External_Face.CellNeighborIndex(i) = str2double( Temp{5} ); % An external face has only one connection (to only one cell).
                External_Face.CellNeighborVec(i,:) = str2double( strsplit( Temp{6} , ';' ) );
            end
            fprintf('\n');
            CornerPointGridData.External_Face = External_Face;
            
            % Reading cell information line by line
            temp = strfind(obj.CornerPointGridMatrix, 'CELL_GEOMETRY');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "CELL_GEOMETRY" is missing. Please check the CornerPointGrid input file!\n');
            end
            Cell.NW_Top_Corner = zeros(CornerPointGridData.N_ActiveCells , 3);
            Cell.NE_Top_Corner = zeros(CornerPointGridData.N_ActiveCells , 3);
            Cell.SW_Top_Corner = zeros(CornerPointGridData.N_ActiveCells , 3);
            Cell.SE_Top_Corner = zeros(CornerPointGridData.N_ActiveCells , 3);
            Cell.NW_Bot_Corner = zeros(CornerPointGridData.N_ActiveCells , 3);
            Cell.NE_Bot_Corner = zeros(CornerPointGridData.N_ActiveCells , 3);
            Cell.SW_Bot_Corner = zeros(CornerPointGridData.N_ActiveCells , 3);
            Cell.SE_Bot_Corner = zeros(CornerPointGridData.N_ActiveCells , 3);
            Cell.Centroid = zeros(CornerPointGridData.N_ActiveCells , 3);
            Cell.Volume   = zeros(CornerPointGridData.N_ActiveCells , 1);
            Cell.N_Neighbors = zeros(CornerPointGridData.N_ActiveCells , 1);
            Cell.Index_Neighbors = cell(CornerPointGridData.N_ActiveCells , 1);
            
            fprintf('---> Cell ');
            for i = 1 : CornerPointGridData.N_ActiveCells
                if (i>1),  fprintf(repmat('\b', 1, 13));  end
                fprintf('%06d/%06d',i,CornerPointGridData.N_ActiveCells);
                Temp = strsplit(obj.CornerPointGridMatrix{index+2-1+i},{' , '});
                Cell.NW_Top_Corner(i,:) = str2double( strsplit( Temp{2} , ';' ) );
                Cell.NE_Top_Corner(i,:) = str2double( strsplit( Temp{3} , ';' ) );
                Cell.SW_Top_Corner(i,:) = str2double( strsplit( Temp{4} , ';' ) );
                Cell.SE_Top_Corner(i,:) = str2double( strsplit( Temp{5} , ';' ) );
                Cell.NW_Bot_Corner(i,:) = str2double( strsplit( Temp{6} , ';' ) );
                Cell.NE_Bot_Corner(i,:) = str2double( strsplit( Temp{7} , ';' ) );
                Cell.SW_Bot_Corner(i,:) = str2double( strsplit( Temp{8} , ';' ) );
                Cell.SE_Bot_Corner(i,:) = str2double( strsplit( Temp{9} , ';' ) );
                Cell.Centroid(i,:) = str2double( strsplit( Temp{10} , ';' ) );
                Cell.Volume(i) = str2double( Temp{11} );
                % Obtaining the cell neighbors
                faceIndex = find(Internal_Face.CellNeighbor1Index==i);
                neighborIndex = Internal_Face.CellNeighbor2Index(faceIndex);
                faceIndex = find(Internal_Face.CellNeighbor2Index==i);
                neighborIndex = [neighborIndex; Internal_Face.CellNeighbor1Index(faceIndex)];
                Cell.N_Neighbors(i) = length(neighborIndex);
                Cell.Index_Neighbors{i} = sort(neighborIndex);
            end
            fprintf('\n');
            CornerPointGridData.Cell = Cell;
        end
        function CornerPointGridRockPropertiesData = ReadCornerPointGridRockPropertiesData(obj)
            fprintf('Reading the CornerPointGridRockProperties input file:\n');
            
            temp = strfind(obj.CornerPointGridRockPropertiesMatrix, 'RESERVOIR_GRID_NX');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "RESERVOIR_GRID_NX" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridRockPropertiesData.Nx = str2double( obj.CornerPointGridRockPropertiesMatrix{index+1} );
            
            temp = strfind(obj.CornerPointGridRockPropertiesMatrix, 'RESERVOIR_GRID_NY');
            if isempty(index)
                error('The keyword "RESERVOIR_GRID_NY" is missing. Please check the CornerPointGrid input file!\n');
            end
            index = find(~cellfun('isempty', temp));
            CornerPointGridRockPropertiesData.Ny = str2double( obj.CornerPointGridRockPropertiesMatrix{index+1} );
                
            temp = strfind(obj.CornerPointGridRockPropertiesMatrix, 'RESERVOIR_GRID_NZ');
            if isempty(index)
                error('The keyword "RESERVOIR_GRID_NZ" is missing. Please check the CornerPointGrid input file!\n');
            end
            index = find(~cellfun('isempty', temp));
            CornerPointGridRockPropertiesData.Nz = str2double( obj.CornerPointGridRockPropertiesMatrix{index+1} );
                
            temp = strfind(obj.CornerPointGridRockPropertiesMatrix, 'ACTIVE_CELLS');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "ACTIVE_CELLS" is missing. Please check the CornerPointGrid input file!\n');
            end
            CornerPointGridRockPropertiesData.N_ActiveCells = str2double( obj.CornerPointGridRockPropertiesMatrix{index+1} );
                
            % Reading Permeability values information line by line(Kx,Ky,Kz)
            temp = strfind(obj.CornerPointGridRockPropertiesMatrix, 'PERMEABILITY_XYZ');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "PERMEABILITY_XYZ" is missing. Please check the CornerPointGridRockProperties input file!\n');
            end
            
            Permeability.X_Y_Z = zeros(CornerPointGridRockPropertiesData.N_ActiveCells, 3);
                
            fprintf('---> Permeability X Y Z ');
            for i = 1 : CornerPointGridRockPropertiesData.N_ActiveCells
%                 if (i>1),  fprintf(repmat('\b', 1, 13));
%                 end
%                 fprintf('%06d/%06d',i,ReservoirProperties.CornerPointGridRockPropertiesData.N_ActiveCells);
                Temp = strsplit(obj.CornerPointGridRockPropertiesMatrix{index+2-1+i},{' , '});
                Permeability.X_Y_Z(i,:) = str2double( strsplit( Temp{2} , ';' ) );
            end
            fprintf('\n');
                       
            temp = strfind(obj.CornerPointGridRockPropertiesMatrix, 'PERMEABILITY_TENSOR');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "PERMEABILITY_TENSOR" is missing. Please check the CornerPointGridRockProperties input file!\n');
            end
                
            Permeability.Tensor = zeros(CornerPointGridRockPropertiesData.N_ActiveCells, 9);
                
            fprintf('---> Permeability Tensor ');
            for i = 1 : CornerPointGridRockPropertiesData.N_ActiveCells
%                 if (i>1),  fprintf(repmat('\b', 1, 13));
%                 end
%                 fprintf('%06d/%06d',i,ReservoirProperties.CornerPointGridRockPropertiesData.N_ActiveCells);
                Temp = strsplit(obj.CornerPointGridRockPropertiesMatrix{index+2-1+i},{' , '});
                Permeability.Tensor(i,1:3) = str2double( strsplit( Temp{2} , ';' ) );
                Permeability.Tensor(i,4:6) = str2double( strsplit( Temp{3} , ';' ) );
                Permeability.Tensor(i,7:9) = str2double( strsplit( Temp{4} , ';' ) );
            end
            fprintf('\n');
            
            CornerPointGridRockPropertiesData.Permeability = Permeability;
            
            % Reading Porosity values data line by line
            temp = strfind(obj.CornerPointGridRockPropertiesMatrix, 'POROSITY');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "POROSITY" is missing. Please check the CornerPointGridRockProperties input file!\n');
            end
                                   
            fprintf('---> Porosity ');
            Porosity = zeros(CornerPointGridRockPropertiesData.N_ActiveCells, 1);

            for i = 1 : CornerPointGridRockPropertiesData.N_ActiveCells
%                 if (i>1),  fprintf(repmat('\b', 1, 13));
%                 end
%                 fprintf('%06d/%06d',i,ReservoirProperties.CornerPointGridRockPropertiesData.N_ActiveCells);
                Temp = strsplit(obj.CornerPointGridRockPropertiesMatrix{index+2-1+i},{' , '});
                Porosity(i) = str2double(Temp{2});
            end
            fprintf('\n');
           
            CornerPointGridRockPropertiesData.Porosity = Porosity;
                        
        end
        function FluidProperties = ReadFluidProperties(obj)
            %%%%%%%%%%%%%FLUID PROPERTIES%%%%%%%%%%%%%%%%
            % 1. Fluid model
            temp = strfind(obj.InputMatrix, 'FLUID MODEL');
            index = find(~cellfun('isempty', temp));
            FluidProperties.FluidModel = char(obj.InputMatrix{index+1});
            FluidProperties.NofPhases = str2double(obj.InputMatrix{index+3});
            FluidProperties.NofComponents = str2double(obj.InputMatrix{index+5});
            if FluidProperties.FluidModel == "Geothermal_2T"
                temp = strfind(obj.InputMatrix, 'AVERAGED_TEMPERATURE');
                T_ave = find(~cellfun('isempty', temp));
                if isempty(T_ave)
                    FluidProperties.AveragedTemperature = "Off";
                else
                    FluidProperties.AveragedTemperature = "On";
                end
            end
            % 2. Density
            temp = strfind(obj.InputMatrix, 'DENSITY');
            index_density = find(~cellfun('isempty', temp));
            % 3. Viscosity
            temp = strfind(obj.InputMatrix, 'VISCOSITY');
            index_viscosity = find(~cellfun('isempty', temp));
            % 5. Liquid Compressibility
            temp = strfind(obj.InputMatrix, 'COMPRESSIBILITY');
            index_comp = find(~cellfun('isempty', temp));
            % 6. Conductivity
            temp = strfind(obj.InputMatrix, 'CONDUCTIVITY');
            index_conduc = find(~cellfun('isempty', temp));
            % 7. Specific Heat
            temp = strfind(obj.InputMatrix, 'SPECIFIC HEAT');
            index_spec_heat = find(~cellfun('isempty', temp));
            
            for i=1:FluidProperties.NofPhases
                FluidProperties.Density(i) = str2double(char(obj.InputMatrix{index_density+(i*2)}));
                FluidProperties.mu(i) = str2double(char(obj.InputMatrix{index_viscosity+(i*2)}));
                FluidProperties.Comp(i) = str2double(char(obj.InputMatrix{index_comp+(i*2)}));
                FluidProperties.FluidConductivity(i) = str2double(char(obj.InputMatrix(index_conduc+(i*2))));
                FluidProperties.SpecificHeat(i) = str2double(char(obj.InputMatrix(index_spec_heat+(i*2))));
            end
         
            % 8. Relative permeability
            temp = strfind(obj.InputMatrix, 'RELPERM');
            index = find(~cellfun('isempty', temp));
            FluidProperties.RelPerm.name = char(obj.InputMatrix(index + 1));
            FluidProperties.RelPerm.s_irr(1) = str2double(obj.InputMatrix(index + 3));
            FluidProperties.RelPerm.s_irr(2) = str2double(obj.InputMatrix(index + 5));
            FluidProperties.RelPerm.Kre(1) = str2double(obj.InputMatrix(index + 7));
            FluidProperties.RelPerm.Kre(2) = str2double(obj.InputMatrix(index + 9));
            FluidProperties.RelPerm.n(1) = str2double(obj.InputMatrix(index + 11));
            FluidProperties.RelPerm.n(2) = str2double(obj.InputMatrix(index + 13));

            % 9. Capillarity
            temp = strfind(obj.InputMatrix, 'CAPILLARITY');
            index = find(~cellfun('isempty', temp));
            FluidProperties.Capillarity.name = char(obj.InputMatrix(index + 1));
            if ~isempty(FluidProperties.Capillarity.name) 
                FluidProperties.Capillarity.wetting = str2double(obj.InputMatrix(index + 2));
            end
           
            % 10. Component properties
            temp = strfind(obj.InputMatrix, 'COMPONENT PROPERTIES');
            index = find(~cellfun('isempty', temp));
            for i = 1:FluidProperties.NofComponents
                FluidProperties.ComponentProperties = str2double(strsplit(char(obj.InputMatrix(index + i * 2))));
            end
            
            % 11. Gravity
            temp = strfind(obj.InputMatrix, 'GRAVITY');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                FluidProperties.Gravity = 'OFF';
            else
                FluidProperties.Gravity = char(obj.InputMatrix(index + 1));
            end
            
        end
        function WellsInfo = ReadWellsInfo(obj, SimulationInput)
            %%%%%%%%%%%%%WELLS%%%%%%%%%%%%%%%%
            temp = regexp(obj.InputMatrix, 'WELL_START', 'match');
            well_start = find(~cellfun('isempty', temp));
            temp = regexp(obj.InputMatrix, 'WELL_END', 'match');
            well_end = find(~cellfun('isempty', temp));
            WellsInfo.NofWell = length(well_start);
            temp  = regexp(obj.InputMatrix, 'INJ', 'match'); 
            inj = find(~cellfun('isempty', temp));
            WellsInfo.NofInj = length(inj);
            temp = regexp(obj.InputMatrix, 'PROD', 'match');
            prod = find(~cellfun('isempty', temp));
            WellsInfo.NofProd = length(prod);
            
            inj=0; prod=0;
            for w = 1 : WellsInfo.NofWell
                WellInputMatrix = obj.InputMatrix(well_start(w)+1:well_end(w)-1);
                temp = regexp(WellInputMatrix, 'TYPE', 'match');
                type = find(~cellfun('isempty', temp));
                temp = regexp(WellInputMatrix, 'COORDINATE', 'match');
                coordinate = find(~cellfun('isempty', temp));
                temp = regexp(WellInputMatrix, 'CONSTRAINT', 'match');
                constraint = find(~cellfun('isempty', temp));
                temp = regexp(WellInputMatrix, 'FORMULA', 'match');
                formula = find(~cellfun('isempty', temp));
                temp = regexp(WellInputMatrix, 'TEMPERATURE', 'match');
                temperature = find(~cellfun('isempty', temp));
                
                % Reading the coordinates of well trajectory
                if strcmp(WellInputMatrix(coordinate+1),'IJK')
                    ijk_1 = strsplit(WellInputMatrix{coordinate+2}, {'	',' ',',','[',']'} );
                    ijk_1 = ijk_1(2:end-1);
                    ijk_1 = strrep(ijk_1,'NX',num2str(SimulationInput.ReservoirProperties.Grid.N(1)));
                    ijk_1 = strrep(ijk_1,'NY',num2str(SimulationInput.ReservoirProperties.Grid.N(2)));
                    ijk_1 = strrep(ijk_1,'NZ',num2str(SimulationInput.ReservoirProperties.Grid.N(3)));
                    ijk_1 = [str2num(ijk_1{1}), str2num(ijk_1{2}), str2num(ijk_1{3})];
                    ijk_2 = strsplit(WellInputMatrix{coordinate+3}, {'	',' ',',','[',']'} );
                    ijk_2 = ijk_2(2:end-1);
                    ijk_2 = strrep(ijk_2,'NX',num2str(SimulationInput.ReservoirProperties.Grid.N(1)));
                    ijk_2 = strrep(ijk_2,'NY',num2str(SimulationInput.ReservoirProperties.Grid.N(2)));
                    ijk_2 = strrep(ijk_2,'NZ',num2str(SimulationInput.ReservoirProperties.Grid.N(3)));
                    ijk_2 = [str2num(ijk_2{1}), str2num(ijk_2{2}), str2num(ijk_2{3})];
                    Well.Coord = [ijk_1;ijk_2];
                    if any(mod(Well.Coord(:),1) ~= 0)
                        error('In well #%d, the ijk coordinates result in non-integer cell indeces. Check the input file!\n', w);
                    end
                elseif strcmp(WellInputMatrix(coordinate+1),'XYZ')
                    Well.Coord = zeros(constraint-coordinate-2,3);
                    for p = 1 : size(Well.Coord,1)
                        xyz = strsplit(WellInputMatrix{coordinate+2+p-1}, {'	',' ',',','[',']'} );
                        xyz = xyz(2:end-1);
                        xyz = strrep(xyz,'LX',num2str(SimulationInput.ReservoirProperties.size(1)));
                        xyz = strrep(xyz,'LY',num2str(SimulationInput.ReservoirProperties.size(2)));
                        xyz = strrep(xyz,'LZ',num2str(SimulationInput.ReservoirProperties.size(3)));
                        Well.Coord(p,:) = [str2num(xyz{1}), str2num(xyz{2}), str2num(xyz{3})];
                    end
                else
                    error('In well #%d, the coordination keyword "IJK" or "XYZ" is missing. Check the input file!\n', w);
                end
                
                % Reading the contraint of well
                Well.Constraint.name = char(WellInputMatrix(constraint+1));
                Well.Constraint.value = str2double(WellInputMatrix(constraint+2));
                
                % Reading the formula type of well
                Well.PI.type = char(WellInputMatrix(formula+1));
                Well.PI.value = str2double(WellInputMatrix(formula+2));
                
                % Reading the temperature of well (only for injection)
                switch SimulationInput.FluidProperties.FluidModel
                    case{'Geothermal_SinlgePhase','Geothermal_MultiPhase'}
                        Well.Temperature = str2double(WellInputMatrix(temperature+1));
                end
                if strcmp(WellInputMatrix(type+1),'INJ')
                    inj = inj+1;
                    WellsInfo.Inj(inj) = Well;
                elseif strcmp(WellInputMatrix(type+1),'PROD')
                    prod = prod+1;
                    WellsInfo.Prod(prod) = Well;
                else
                    error('In well #%d, the type of the well should either be "INJ" or "PROD". Check the input file!\n', w);
                end
            end
        end
        function SimulatorSettings = ReadSimulatorSettings(obj, SimulationInput)
            % 1. Max number of time-steps
            temp = strfind(obj.SettingsMatrix, 'TIMESTEPS');
            xv = find(~cellfun('isempty', temp));
            SimulatorSettings.MaxNumTimeSteps = str2double(obj.SettingsMatrix(xv+1));
            if  strcmp(obj.SettingsMatrix(xv+2),'Fix')
                SimulatorSettings.FixedStep = 1;
            else
                SimulatorSettings.FixedStep = 0;
            end
                % 2. Min and max time-step sizes
            temp = strfind(obj.SettingsMatrix, 'MINMAXDT');
            xv = find(~cellfun('isempty', temp));
            if ~isempty(xv)
                SimulatorSettings.MinMaxdt = str2double(strsplit(char(obj.SettingsMatrix(xv+1))));
            else
                % default value
                SimulatorSettings.MinMaxdt = [0.1, 30];
            end
            temp = strfind(obj.SettingsMatrix, '');
            % 3. Number of reports --> report is written every
            % T_total/n_reports s; if 0 report at every dt
            temp = strfind(obj.SettingsMatrix, 'REPORTS');
            xv = find(~cellfun('isempty', temp));
            SimulatorSettings.reports = str2double(obj.SettingsMatrix(xv+1));
            % 4. Choose coupling
            temp = strfind(obj.SettingsMatrix, 'FIM'); 
            index = find(~cellfun('isempty', temp));
            if  ~isempty(index)
                SimulatorSettings.CouplingType = 'FIM';
                SimulatorSettings.MaxIterations = str2double(obj.SettingsMatrix(index + 1)); 
                SimulatorSettings.ResidualTolerances = str2double(strsplit(regexprep(obj.SettingsMatrix{index + 2},' ' ,''),','));
                SimulatorSettings.SolutionTolerances = str2double(strsplit(regexprep(obj.SettingsMatrix{index + 3},' ' ,''),','));
                SimulatorSettings.cfl = str2double(obj.SettingsMatrix(index + 4));
            else
                SimulatorSettings.CouplingType ='Sequential';
                temp = strfind(obj.SettingsMatrix, 'SEQUENTIAL'); 
                index = find(~cellfun('isempty', temp));
                SimulatorSettings.MaxIterations = str2double(obj.SettingsMatrix(index + 1));
                SimulatorSettings.Tolerance     = str2double(obj.SettingsMatrix(index + 2));
                SimulatorSettings.cfl           = str2double(obj.SettingsMatrix(index + 3));
                SimulatorSettings.TransportSolver.Type = 'EXPLICIT';
                SimulatorSettings.TransportSolver.Tol  = 1e-3;
                SimulatorSettings.TransportSolver.MaxIter  = 100;
                temp = strfind(obj.SettingsMatrix, 'IMPSAT');
                index = find(~cellfun('isempty', temp));
                if ~isempty(index)
                    SimulatorSettings.TransportSolver.Type = 'IMPSAT';
                    SimulatorSettings.TransportSolver.Tol  = str2double(obj.SettingsMatrix(index + 2));
                    SimulatorSettings.TransportSolver.MaxIter  = str2double(obj.SettingsMatrix(index + 3));
                end
            end
            
            % 5. Choose formulation
            temp = strfind(obj.SettingsMatrix, 'FORMULATION');
            xv = find(~cellfun('isempty', temp));
            if isempty(xv)
                switch (SimulationInput.FluidProperties.FluidModel)
                    case('SinglePhase')
                        SimulatorSettings.Formulation = 'Immiscible';
                    case('Immiscible')
                        SimulatorSettings.Formulation = 'Immiscible';
                    case("Geothermal_1T")
                        SimulatorSettings.Formulation = "Geothermal_1T";
                    case("Geothermal_2T")
                        SimulatorSettings.Formulation = "Geothermal_2T";
                    otherwise
                        SimulatorSettings.Formulation = 'Molar';
                end
            else
                SimulatorSettings.Formulation =  char(obj.SettingsMatrix(xv+1));
            end
            if strcmp(SimulatorSettings.CouplingType, 'SEQUENTIAL')
                SimulatorSettings.Formulation = 'Immiscible';
            end
            
            % 6. Linear solver type 
            temp = strfind(obj.SettingsMatrix, 'LINEARSOLVER'); 
            xv = find(~cellfun('isempty', temp));
            SimulatorSettings.LinearSolver =  char(obj.SettingsMatrix(xv+1));
            
            % Default value is fine-scale 
            SimulatorSettings.DiscretizationModel = 'FS';
            
            %% ADM settings
            temp = strfind(obj.SettingsMatrix, 'ADM');
            adm = find(~cellfun('isempty', temp));
            if str2double(obj.SettingsMatrix(adm + 1)) == 1
                SimulatorSettings.DiscretizationModel = 'ADM';
                % ADM settings for reservoir
                temp = strfind(obj.SettingsMatrix, 'LEVELS');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.ADMSettings.maxlevel = zeros(1+SimulationInput.FracturesProperties.NrOfFrac, 1);
                SimulatorSettings.ADMSettings.maxLevel(1) = str2double(obj.SettingsMatrix(x+1));
                SimulatorSettings.ADMSettings.Coarsening = zeros( 1+SimulationInput.FracturesProperties.NrOfFrac, 3, SimulatorSettings.ADMSettings.maxLevel(1) );
                temp = strfind(obj.SettingsMatrix, 'COARSENING_RATIOS');
                x = find(~cellfun('isempty', temp));
                cx = str2double(obj.SettingsMatrix(x+1));
                cy = str2double(obj.SettingsMatrix(x+2));
                cz = str2double(obj.SettingsMatrix(x+3));
                temp = strfind(obj.SettingsMatrix, 'COARSENING_CRITERION');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.ADMSettings.GridSelCriterion = char(obj.SettingsMatrix(x+1));
                temp = strfind(obj.SettingsMatrix, 'VARIABLE');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.ADMSettings.key = char(obj.SettingsMatrix(x+1));
                temp = strfind(obj.SettingsMatrix, 'TOLERANCE');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.ADMSettings.tol = str2double(obj.SettingsMatrix(x+1));
                for L = 1:SimulatorSettings.ADMSettings.maxLevel(1)
                    SimulatorSettings.ADMSettings.Coarsening(1,:,L) = [cx, cy, cz].^L; %Coarsening Factors: Cx1, Cy1; Cx2, Cy2; ...; Cxn, Cyn;
                end
                temp = strfind(obj.SettingsMatrix, 'PRESSURE_INTERPOLATOR');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.ADMSettings.PInterpolator = char(obj.SettingsMatrix(x+1));
                temp = strfind(obj.SettingsMatrix, 'HYPERBOLIC_INTERPOLATOR');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.ADMSettings.HInterpolator = char(obj.SettingsMatrix(x+1));
                temp = strfind(obj.SettingsMatrix, 'ROCK_TEMPERATURE_INTERPOLATOR');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.ADMSettings.TrInterpolator = char(obj.SettingsMatrix(x+1));

                if isempty(SimulatorSettings.ADMSettings.maxLevel(1)) || isempty(SimulatorSettings.ADMSettings.Coarsening(1,:,:)) || isempty(SimulatorSettings.ADMSettings.key) || isempty(SimulatorSettings.ADMSettings.tol) || isempty(SimulatorSettings.ADMSettings.PInterpolator)
                    error('DARSIM2 ERROR: Missing ADM settings! Povide LEVELS, COARSENING_CRITERION, COARSENING_RATIOS, TOLERANCE, PRESSURE_INTERPOLATOR');
                end
                
                temp = strfind(obj.SettingsMatrix, 'DLGR');
                x = find(~cellfun('isempty', temp));
                if isempty(temp)
                    SimulatorSettings.ADMSettings.DLGR = 0;
                else
                    SimulatorSettings.ADMSettings.DLGR = str2double(obj.SettingsMatrix(x+1));
                end
                
                % Coupled or de-coupled basis function
                temp = strfind(obj.SettingsMatrix, 'PRESSURE_INTERPOLATOR');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.ADMSettings.BFtype = char(obj.SettingsMatrix(x+2));
                if ~strcmp(SimulatorSettings.ADMSettings.BFtype,'COUPLED') && ~strcmp(SimulatorSettings.ADMSettings.BFtype,'DECOUPLED')
                    SimulatorSettings.ADMSettings.BFtype= 'COUPLED';
                end
                
                % Maximum heterogeneity contrast for basis functions
                temp = strfind(obj.SettingsMatrix, 'BASIS_FUNCTION_MAX_CONTRAST');
                Index = find(~cellfun('isempty', temp));
                if isempty(Index)
                    SimulatorSettings.ADMSettings.BF_MaxContrast = 1e2;
                else
                    SimulatorSettings.ADMSettings.BF_MaxContrast = str2double(obj.SettingsMatrix(Index+1));
                end
                
                % Maximum heterogeneity contrast of pEDFM alpha factors for basis functions
                temp = strfind(obj.SettingsMatrix, 'pEDFM_MAX_CONTRAST');
                Index = find(~cellfun('isempty', temp));
                if isempty(Index)
                    SimulatorSettings.ADMSettings.pEDFM_MaxContrast = 1e1;
                else
                    SimulatorSettings.ADMSettings.pEDFM_MaxContrast = str2double(obj.SettingsMatrix(Index+1));
                end
                
                % ADM settings in the fractures
                if SimulationInput.FracturesProperties.Fractured
                    NrOfFrac = SimulationInput.FracturesProperties.NrOfFrac;
                    temp = strfind(obj.FractureMatrix, 'PROPERTIES');
                    frac_index = find(~cellfun('isempty', temp));
                    for f = 1 : NrOfFrac
                        frac_info_split = strsplit(obj.FractureMatrix{frac_index(f)},' ');
                        ADM_temp = regexprep(frac_info_split{9},' ' ,'');
                        ADM_temp = strsplit(ADM_temp, { '[' , ',' , ']' });
                        ADM_temp = [ str2double(ADM_temp(2)) , str2double(ADM_temp(3)) , str2double(ADM_temp(4)) , str2double(ADM_temp(5)) ];
                        if ADM_temp(1)
                            SimulatorSettings.ADMSettings.maxLevel(1+f) = ADM_temp(2);
                        else
                            SimulatorSettings.ADMSettings.maxLevel(1+f) = 0;
                        end
                        for L = 1:SimulatorSettings.ADMSettings.maxLevel(1)
                            if L <= SimulatorSettings.ADMSettings.maxLevel(1+f)
                                SimulatorSettings.ADMSettings.Coarsening(1+f,:,L) = [ADM_temp(3), ADM_temp(4), 1].^L;
                            else
                                SimulatorSettings.ADMSettings.Coarsening(1+f,:,L) = [ADM_temp(3), ADM_temp(4), 1].^SimulatorSettings.ADMSettings.maxLevel(1+f);
                            end
                        end
                    end
                end
            end
            
           %% Multilevel Multiscale settings
            temp = strfind(obj.SettingsMatrix, 'MMs');
            mms = find(~cellfun('isempty', temp));
            if str2double(obj.SettingsMatrix(mms + 1)) == 1
                SimulatorSettings.DiscretizationModel = 'MMs';
                temp = strfind(obj.SettingsMatrix, 'LEVELS');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.MMsSettings.maxlevel = zeros(1+SimulationInput.FracturesProperties.NrOfFrac, 1);
                SimulatorSettings.MMsSettings.maxLevel(1) = str2double(obj.SettingsMatrix(x+1));
                SimulatorSettings.MMsSettings.Coarsening = zeros( 1+SimulationInput.FracturesProperties.NrOfFrac, 3, SimulatorSettings.MMsSettings.maxLevel(1) );
                temp = strfind(obj.SettingsMatrix, 'COARSENING_RATIOS');
                x = find(~cellfun('isempty', temp));
                cx = str2double(obj.SettingsMatrix(x+1));
                cy = str2double(obj.SettingsMatrix(x+2));
                cz = str2double(obj.SettingsMatrix(x+3));
                for L = 1:SimulatorSettings.MMsSettings.maxLevel(1)
                    SimulatorSettings.MMsSettings.Coarsening(1,:,L) = [cx, cy, cz].^L; %Coarsening Factors: Cx1, Cy1; Cx2, Cy2; ...; Cxn, Cyn;
                end
                
                % If you write any of these keywords in the input file the
                % options will be active (thus default values are 0)
                temp = sum(contains(obj.SettingsMatrix, 'MSFE'));
                if temp
                    SimulatorSettings.MMsSettings.MSFE = 1;
                else
                    SimulatorSettings.MMsSettings.MSFE = 0;
                end
                
                temp = sum(contains(obj.SettingsMatrix, 'CORRECTION FUNCTIONS'));
                if temp
                    SimulatorSettings.MMsSettings.CorrectionFunctions = 1;
                else
                    SimulatorSettings.MMsSettings.CorrectionFunctions = 0;
                end

                % Coupled or de-coupled basis function
                temp = strfind(obj.SettingsMatrix, 'PRESSURE_INTERPOLATOR');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.MMsSettings.BFtype = char(obj.SettingsMatrix(x+2));
                if ~strcmp(SimulatorSettings.MMsSettings.BFtype,'COUPLED') && ~strcmp(SimulatorSettings.MMsSettings.BFtype,'DECOUPLED')
                    SimulatorSettings.MMsSettings.BFtype= 'COUPLED';
                end
                
                % Maximum heterogeneity contrast for basis functions
                temp = strfind(obj.SettingsMatrix, 'BASIS_FUNCTION_MAX_CONTRAST');
                Index = find(~cellfun('isempty', temp));
                if isempty(Index)
                    SimulatorSettings.MMsSettings.BF_MaxContrast = 1e2;
                else
                    SimulatorSettings.MMsSettings.BF_MaxContrast = str2double(obj.SettingsMatrix(Index+1));
                end
                
                % Maximum heterogeneity contrast of pEDFM alpha factors for basis functions
                temp = strfind(obj.SettingsMatrix, 'pEDFM_MAX_CONTRAST');
                Index = find(~cellfun('isempty', temp));
                if isempty(Index)
                    SimulatorSettings.MMsSettings.pEDFM_MaxContrast = 1e1;
                else
                    SimulatorSettings.MMsSettings.pEDFM_MaxContrast = str2double(obj.SettingsMatrix(Index+1));
                end
                
                % MMs settings in fractures
                if SimulationInput.FracturesProperties.Fractured
                    NrOfFrac = SimulationInput.FracturesProperties.NrOfFrac;
                    temp = strfind(obj.FractureMatrix, 'PROPERTIES');
                    frac_index = find(~cellfun('isempty', temp));
                    for f = 1 : NrOfFrac
                        frac_info_split = strsplit(obj.FractureMatrix{frac_index(f)},' ');
                        MMs_temp = regexprep(frac_info_split{9},' ' ,'');
                        MMs_temp = strsplit(MMs_temp, { '[' , ',' , ']' });
                        MMs_temp = [ str2double(MMs_temp(2)) , str2double(MMs_temp(3)) , str2double(MMs_temp(4)) , str2double(MMs_temp(5)) ];
                        if MMs_temp(1)
                            SimulatorSettings.MMsSettings.maxLevel(1+f) = MMs_temp(2);
                        else
                            SimulatorSettings.MMsSettings.maxLevel(1+f) = 0;
                        end
                        for L = 1:SimulatorSettings.MMsSettings.maxLevel(1)
                            if L <= SimulatorSettings.MMsSettings.maxLevel(1+f)
                                SimulatorSettings.MMsSettings.Coarsening(1+f,:,L) = [MMs_temp(3), MMs_temp(4), 1].^L;
                            else
                                SimulatorSettings.MMsSettings.Coarsening(1+f,:,L) = [MMs_temp(3), MMs_temp(4), 1].^SimulatorSettings.MMsSettings.maxLevel(1+f);
                            end
                        end
                    end
                end
            end
            
            %% %%% LTS
            % temp = sum(contains(obj.SettingsMatrix, 'LTS'));
			lts = strfind(obj.SettingsMatrix, 'LTS');
            index = find(~cellfun('isempty', lts));
            if str2double(obj.SettingsMatrix(index + 1)) == 1
                SimulatorSettings.LTS = 1;
                temp = strfind(obj.SettingsMatrix, 'REF_CRITERION');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.LTSCriterion = char(obj.SettingsMatrix(x+1));
                temp = strfind(obj.SettingsMatrix, 'TIME_TOL');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.LTStol = str2double(obj.SettingsMatrix(x+1));
                temp2 = strfind(obj.SettingsMatrix, 'PLOT_INTERNAL_STEPS');
                x2 = find(~cellfun('isempty', temp2));
                SimulatorSettings.LTSPlot = str2double(obj.SettingsMatrix(x2 + 1));
            else
                SimulatorSettings.LTS = 0;
                SimulatorSettings.LTSPlot = 0;
            end
            

            
            %% %%% Stop criterion
            SimulatorSettings.StopCriterion = 'MAX TIME'; % Decide up to when you want to run the simulation
            
            %% %%%%%%%%%%%OPTIONS%%%%%%%%%%%%%%%%
            temp = strfind(obj.SettingsMatrix, 'OUTPUT'); 
            xv = find(~cellfun('isempty', temp));
            SimulatorSettings.plotting = char(obj.SettingsMatrix(xv+1)); %Matlab or VTK
        end
    end
end