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
        function obj = reader_darsim2(dir, file)
            obj@reader(dir, file);          
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
            temp = strfind(obj.InputMatrix, 'FRACTURED'); % Check the main input file and look for FRACTURED
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                Boolean = obj.InputMatrix{index+1};
                if strcmp(Boolean,'ON')
                    FractureFile = strcat(obj.Directory,'/',obj.InputMatrix{index+2});
                    if ~isfile(FractureFile)
                        error('DARSim Error: The fracture input file "%s" does not exist!\n', FractureFile);
                    end
                    fileID = fopen(FractureFile, 'r');
                    fprintf('Reading fracture file ...');
                    matrix = textscan(fileID, '%s', 'Delimiter', '\n');
                    obj.FractureMatrix = matrix{1};
                    fclose(fileID);
                    fprintf('--> Successful.\n\n');
                end
            end

            [Builder.SimulationInput, Builder.SimulatorSettings] = ReadInformation(obj);
        end
        function [SimulationInput, SimulatorSettings] = ReadInformation(obj)
            %% 1. INPUT FILE
            % 1.0. Passinbg the main input directory to simulation_builder
            SimulationInput.Directory = obj.Directory;
            
            % 1.1. Name of the Problem 
            temp = strfind(obj.InputMatrix, 'TITLE'); % Search a specific string and find all rows containing matches
            SimulationInput.ProblemName = char(obj.InputMatrix(find(~cellfun('isempty', temp)) + 1));
			
            % 1.2. Total time
            temp = strfind(obj.InputMatrix, 'TOTAL_TIME');
            xv = find(~cellfun('isempty', temp));
            if isempty(xv)
                error('The keyword "TOTAL_TIME" is missing. Please check the input file!\n');
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
			
            % 1.3. Reservoir properties
            SimulationInput.ReservoirProperties = obj.ReadReservoirProperties();
            
            % 1.4. Fluid properties
            SimulationInput.FluidProperties = obj.ReadFluidProperties();
            
            % 1.5. Rock properties
            SimulationInput.RockProperties = obj.ReadRockProperties();

            % 1.6. Initial conditions
            SimulationInput.InitialConditions = obj.ReadInitialConditions();
            SimulationInput.ReservoirProperties.Temperature = SimulationInput.InitialConditions.Temperature;
            
            % 1.7. Wells info
            SimulationInput.WellsInfo = obj.ReadWellsInfo(SimulationInput);
            
            % 1.8. Fractures' properties (if any)
            temp = strfind(obj.InputMatrix, 'FRACTURED'); % Check the main input file and look for FRACTURED
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                SimulationInput.FracturesProperties.isFractured = 0;
                SimulationInput.FracturesProperties.NumOfFrac = 0;
            else
                Boolean = obj.InputMatrix{index+1};
                if strcmp(Boolean,'ON')
                    SimulationInput.FracturesProperties = obj.ReadFracturesProperties();
                else
                    SimulationInput.FracturesProperties.isFractured = 0;
                    SimulationInput.FracturesProperties.NumOfFrac = 0;
                end
            end
            
            %% 2. SIMULATOR'S SETTINGS
            SimulatorSettings = obj.ReadSimulatorSettings(SimulationInput);
        end
        function ReservoirProperties = ReadReservoirProperties(obj)
            %% 1. Discretization Approach
            temp = strfind(obj.InputMatrix, 'GRID_DISCRETIZATION');
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                ReservoirProperties.Discretization = obj.InputMatrix{index+1};
            else
                warning('The keyword "GRID_DISCRETIZATION" is missing. Cartesian discretization is set by default.\n');
                ReservoirProperties.Discretization = 'CartesianGrid';
            end
            
            % Check the discretization approach
            switch ReservoirProperties.Discretization
                case('CornerPointGrid')
                    % Read the cornerpointgrid data input file
                    FileName = obj.InputMatrix{index+2};
                    if isfile(FileName)
                        CornerPointGridFile = FileName;
                    elseif isfile(strcat(obj.Directory, '/', FileName))
                        CornerPointGridFile = strcat(obj.Directory, '/', FileName);
                    else
                        error('The path to CornerPointGrid input file is not valid. Please specify it correctly in the main input file.');
                    end
                    
                    % Check if the rock properties of CornerPointGrid model is set to be used or not
                    temp = strfind(obj.InputMatrix, 'INCLUDE_ROCK_PROPERTIES');
                    index = find(~cellfun('isempty', temp));
                    if ~isempty(index)
                        includeCornerPointGridRockProperties = 1;
                    else
                        includeCornerPointGridRockProperties = 0;
                    end
                    
                    % Check if the data has already been loaded and saved, otherwise read the CornerPointGrid input file
                    if isfile(strcat(obj.Directory, '/','CornerPointGridData.mat'))
                        fprintf('"CornerPointGridData.mat" file already exists. No need to load the CornerPointGrid data input file.\n');
                        fprintf('---> Loading the "CornerPointGridData.mat" file ...');
                        load(strcat(obj.Directory, '/','CornerPointGridData.mat'),'ReservoirProperties');
                        fprintf('Done!\n');
                    else
                        EclipseReader = reader_eclipse(CornerPointGridFile);
                        if contains(CornerPointGridFile,'.grdecl','IgnoreCase',true) && ~contains(CornerPointGridFile,'.txt')
                            % Reading the grdecl file directly
                            ReservoirProperties.CornerPointGridData = EclipseReader.ReadGRDECL(1);
                        elseif contains(CornerPointGridFile,'.txt')
                            % Reading the txt files generated by converting grdecl via DARSim CornerPointGrid data generator   
                            FileName = '';
                            % adding the cornerpointgrid rock properties data input file to ECLIPSE reader (if any)
                            if ~isempty(index)
                                FileName = obj.InputMatrix{index+1};
                            end
                            if isfile(FileName)
                                CornerPointGridRockPropertiesFile = FileName;
                            elseif isfile(strcat(obj.Directory, '/', FileName))
                                CornerPointGridRockPropertiesFile = strcat(obj.Directory, '/', FileName);
                            end
                            ReservoirProperties.CornerPointGridData = EclipseReader.ReadfromTXT(CornerPointGridRockPropertiesFile);
                        else
                            error('The extension of CornerPointGrid input file is not supported (must be txt or grdecl).');
                        end
                        
                        fprintf('---> Saving the "CornerPointGridData.mat" file ...');
                        save(strcat(obj.Directory, '/','CornerPointGridData.mat'),'ReservoirProperties','-v7.3');
                        fprintf('Done!\n');
                    end
                    
                    % Deleting the permeability data if simulation is set to be homogeneous
                    if ~includeCornerPointGridRockProperties
                        ReservoirProperties.CornerPointGridData = rmfield(ReservoirProperties.CornerPointGridData,{'Porosity','Permeability','PermUnit','PermScale'});
                    end
                    
                    ReservoirProperties.Grid.N(1) = ReservoirProperties.CornerPointGridData.Nx;
                    ReservoirProperties.Grid.N(2) = ReservoirProperties.CornerPointGridData.Ny;
                    ReservoirProperties.Grid.N(3) = ReservoirProperties.CornerPointGridData.Nz;
                    ReservoirProperties.Grid.N_ActiveCells = ReservoirProperties.CornerPointGridData.N_ActiveCells;
                    
                    ReservoirProperties.size = [ReservoirProperties.CornerPointGridData.Lx;...
                                                ReservoirProperties.CornerPointGridData.Ly;
                                                ReservoirProperties.CornerPointGridData.Lz];
                    
                case('CartesianGrid')
                    % Assume it is Cartesian
                    % 1. size of the reservoir
                    temp = strfind(obj.InputMatrix, 'DIMENSION'); % Search a specific string and find all rows containing matches
                    index = find(~cellfun('isempty', temp));
                    if isempty(index)
                        error('The keyword "DIMENSION" is missing. Please check the input file!\n');
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
                    ReservoirProperties.Grid.N_ActiveCells = prod(ReservoirProperties.Grid.N);
                otherwise
                    error('The discretization method should either be "CartesianGrid" or "CornerPointGrid". Check the input file!\n');
            end
            
            
            %% 2. Permeability
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
            
            % Read the permeability
            if strcmp(ReservoirProperties.Discretization,'CornerPointGrid') && sum(strcmp(fieldnames(ReservoirProperties.CornerPointGridData), 'Permeability'))
                % Do nothing, beacuse the permeability data already exist in CornerPointGridData.
                ReservoirProperties.Perm = ReservoirProperties.CornerPointGridData.Permeability;
                if ~strcmp(ReservoirProperties.PermUnit, ReservoirProperties.CornerPointGridData.PermUnit)
                    warning('The permeability unit reported in the main input file and the CornerPointGrid rock properties input file do not match.\nThe one from CornerPointGrid rock properties input file (%s) is considered.\n', ReservoirProperties.CornerPointGridData.PermUnit);
                end
                if ~strcmp(ReservoirProperties.PermScale, ReservoirProperties.CornerPointGridData.PermScale)
                    warning('The permeability scale reported in the main input file and the CornerPointGrid rock properties input file do not match.\nThe one from CornerPointGrid rock properties input file (%s) is considered.\n', ReservoirProperties.CornerPointGridData.PermScale);
                end
            else
                % The permeability is read from the main input file.
                ReservoirProperties.Perm = [];
                perm = zeros(3, 1);
                temp = strfind(obj.InputMatrix, 'PERMEABILITY_X');
                perm(1) = find(~cellfun('isempty', temp));
                if isempty(perm(1))
                    error('The keyword "PERMEABILITY_X" is missing. Please check the input file!\n');
                end
                temp = strfind(obj.InputMatrix, 'PERMEABILITY_Y');
                perm(2) = find(~cellfun('isempty', temp));
                if isempty(perm(2))
                    error('The keyword "PERMEABILITY_Y" is missing. Please check the input file!\n');
                end
                temp = strfind(obj.InputMatrix, 'PERMEABILITY_Z');
                perm(3) = find(~cellfun('isempty', temp));
                if isempty(perm(3))
                    error('The keyword "PERMEABILITY_Z" is missing. Please check the input file!\n');
                end
                DimChar = {'x' , 'y' , 'z'};
                for i=1:3
                    if contains(obj.InputMatrix(perm(i) - 1), 'INCLUDE')
                        ReservoirProperties.PermInclude(i) = 1;
                        if isfile( char(obj.InputMatrix(perm(i)+1)) )
                            ReservoirProperties.PermFile{i} = char(obj.InputMatrix(perm(i)+1));
                        elseif isfile( strcat( obj.Directory, '/', char(obj.InputMatrix(perm(i)+1)) ) )
                            ReservoirProperties.PermFile{i} = strcat( obj.Directory, '/', char(obj.InputMatrix(perm(i)+1)) );
                        else
                            error('The path to permeability input file is not valid. Please specify it correctly in the main input file.');
                        end
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
            end
            
            %% 3. Porosity
            if strcmp(ReservoirProperties.Discretization,'CornerPointGrid') && sum(strcmp(fieldnames(ReservoirProperties.CornerPointGridData), 'Porosity'))
                % Do nothing, beacuse the porosity data already exist in CornerPointGridData.
                ReservoirProperties.phi = ReservoirProperties.CornerPointGridData.Porosity;
                ReservoirProperties.PorosityInclude = 0;
            else
                temp = strfind(obj.InputMatrix, 'POROSITY');
                index = find(~cellfun('isempty', temp));
                if isempty(index)
                    error('The keyword "Por" for porosity is missing. Please check the input file!\n');
                end
                if contains(obj.InputMatrix(index - 1), 'INCLUDE')
                    ReservoirProperties.PorosityInclude = 1;
                    ReservoirProperties.PorosityFile = strcat(char(obj.InputMatrix(index+1)));
                else
                    ReservoirProperties.PorosityInclude = 0;
                    ReservoirProperties.phi = str2double(obj.InputMatrix(index + 1));
                end
            end
        end
        function FracturesProperties = ReadFracturesProperties(obj)
            %%%%%%%%%%%%%PROPERTIES OF THE FRACTURES%%%%%%%%%%%%%%%%
            % Marking the simulation as fractured
            FracturesProperties.isFractured = 1;

            % Reading the discretization approach used in fracture generator
            temp = strfind(obj.FractureMatrix, 'DISCRETIZATION');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "DISCRETIZATION" is missing. Please check the Fracture_Output file!\n');
            end
            temp = strsplit(obj.FractureMatrix{index},' ');
            FracturesProperties.DiscretizationType = temp{end};
            
            % Reading the fracture model approach (EDFM or pEDFM)
            temp = strfind(obj.FractureMatrix, 'FRACTURE_MODEL_TYPE');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "FRACTURE_MODEL_TYPE" (for fracture model approach EDFM or pEDFM) is missing. Please check the Fracture_Output file!\n');
            end
            temp = strsplit(obj.FractureMatrix{index},' ');
            FracturesProperties.FractureModelType = temp{end};
            
            % Reading the reservoir size (Lx,Ly,Lz [m3]) from the fracture file
            temp = strfind(obj.FractureMatrix, 'DIMENSION');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "DIMENSION" is missing. Please check the Fracture_Output file!\n');
            end
            temp = strsplit(obj.FractureMatrix{index},' ');
            FracturesProperties.Reservoir.Lx = str2double(temp{2});
            FracturesProperties.Reservoir.Ly = str2double(temp{4});
            FracturesProperties.Reservoir.Lz = str2double(temp{6});
            
            % Reading the reservoir grid (Nx,Ny,Nz) from the fracture file
            temp = strfind(obj.FractureMatrix, 'RESERVOIR_GRID');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "RESERVOIR_GRID" is missing. Please check the Fracture_Output file!\n');
            end
            temp = strsplit(obj.FractureMatrix{index},' ');
            FracturesProperties.Reservoir.Nx = str2double(temp{2});
            FracturesProperties.Reservoir.Ny = str2double(temp{4});
            FracturesProperties.Reservoir.Nz = str2double(temp{6});
            
            % Reading the number of active grids reservoir grids from the fracture file
            temp = strfind(obj.FractureMatrix, 'RESERVOIR_ACTIVE_GRID');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "RESERVOIR_ACTIVE_GRID" is missing. Please check the Fracture_Output file!\n');
            end
            temp = strsplit(obj.FractureMatrix{index},' ');
            FracturesProperties.Reservoir.N = str2double(temp{2});
            
            % Reading the ADM settings of the reservoir from the fracture file
            temp = strfind(obj.FractureMatrix, 'RESERVOIR_ADM');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                warning('The keyword "RESERVOIR_ADM" is missing in the Fracture_Output file.\nThe default values will be set.');
                FracturesProperties.ReservoirADM.isActive = 0;
                FracturesProperties.ReservoirADM.CoarseLevel = 1;
                FracturesProperties.ReservoirADM.CoarseRatio = [ 1 ; 1 ; 1];
            end
            temp = strsplit(obj.FractureMatrix{index},{' ',',','[',']'});
            FracturesProperties.ReservoirADM.isActive = str2double(temp{2});
            FracturesProperties.ReservoirADM.CoarseLevel = str2double(temp{3});
            FracturesProperties.ReservoirADM.CoarseRatio = [ str2double(temp{4}) ; str2double(temp{5}) ; str2double(temp{6}) ];
            
            % Reading the number of fractures
            temp = strfind(obj.FractureMatrix, 'NUM_FRACS');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "NUM_FRACS" is missing. Please check the Fracture_Output file!\n');
            end
            temp = strsplit(obj.FractureMatrix{index},' ');
            FracturesProperties.NumOfFrac = str2double( temp{2} );
            
            % Reading the number of all fractures grids
            temp = strfind(obj.FractureMatrix, 'NUM_TOTAL_FRACS_GRIDS');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                error('The keyword "NUM_TOTAL_FRACS_GRIDS" is missing. Please check the Fracture_Output file!\n');
            end
            temp = strsplit(obj.FractureMatrix{index},' ');
            FracturesProperties.NumOfAllFracGrids = str2double( temp{2} );
        end
        function FluidProperties = ReadFluidProperties(obj)
            % 1. Fluid model
            temp = strfind(obj.InputMatrix, 'FLUID_MODEL');
            index = find(~cellfun('isempty', temp));
            FluidProperties.FluidModel = char(obj.InputMatrix{index+1});
            FluidProperties.NofPhases = str2double(obj.InputMatrix{index+3});
            FluidProperties.NofComponents = str2double(obj.InputMatrix{index+5});
            % 2. Density of the fluid
            temp = strfind(obj.InputMatrix, 'FLUID_DENSITY');
            index_density = find(~cellfun('isempty', temp));
            % 3. Viscosity of the fluid
            temp = strfind(obj.InputMatrix, 'FLUID_VISCOSITY');
            index_viscosity = find(~cellfun('isempty', temp));
            % 4. Liquid Compressibility
            temp = strfind(obj.InputMatrix, 'FLUID_COMPRESSIBILITY');
            index_compressibility = find(~cellfun('isempty', temp));
            % 5. Conductivity
            temp = strfind(obj.InputMatrix, 'FLUID_HEAT_CONDUCTIVITY');
            index_heat_conductivity = find(~cellfun('isempty', temp));
            % 6. Specific Heat
            temp = strfind(obj.InputMatrix, 'FLUID_SPECIFIC_HEAT');
            index_specific_heat = find(~cellfun('isempty', temp));
            
            for i=1:FluidProperties.NofPhases
                FluidProperties.Density(i)           = str2double(char(obj.InputMatrix{index_density           + i}));
                FluidProperties.mu(i)                = str2double(char(obj.InputMatrix{index_viscosity         + i}));
                FluidProperties.Comp(i)              = str2double(char(obj.InputMatrix{index_compressibility   + i}));
                FluidProperties.FluidConductivity(i) = str2double(char(obj.InputMatrix(index_heat_conductivity + i)));
                FluidProperties.SpecificHeat(i)      = str2double(char(obj.InputMatrix(index_specific_heat     + i)));
            end
         
            % 7. Relative permeability
            temp = strfind(obj.InputMatrix, 'RELPERM');
            index = find(~cellfun('isempty', temp));
            FluidProperties.RelPerm.name = char(obj.InputMatrix(index + 1));
            switch FluidProperties.RelPerm.name
                case{'Linear','Quadratic','Corey','Foam'}
                    FluidProperties.RelPerm.s_irr(1) = str2double(obj.InputMatrix(index + 3));
                    FluidProperties.RelPerm.s_irr(2) = str2double(obj.InputMatrix(index + 5));
                    FluidProperties.RelPerm.Kre(1)   = str2double(obj.InputMatrix(index + 7));
                    FluidProperties.RelPerm.Kre(2)   = str2double(obj.InputMatrix(index + 9));
                    FluidProperties.RelPerm.n(1)     = str2double(obj.InputMatrix(index + 11));
                    FluidProperties.RelPerm.n(2)     = str2double(obj.InputMatrix(index + 13));
                case('Table')
                    FluidProperties.RelPerm.FileFormat = char(obj.InputMatrix(index + 3));
                    FluidProperties.RelPerm.TableType  = char(obj.InputMatrix(index + 5));
                    
                    switch FluidProperties.RelPerm.TableType
                        case{'Drainage','Imbibition'}
                            NumOfTable = 1;
                        case('Cyclic')
                            NumOfTable = 2;
                    end
                    for n = 1 : NumOfTable
                        FileName = char(obj.InputMatrix(index + 6+n));
                        FluidProperties.RelPerm.TableFile{n,1} = strcat(obj.Directory,'/',FileName);
                        switch FluidProperties.RelPerm.FileFormat
                            case('MATLAB')
                                % Not implemented yet
                                error('The reading of RelPerm tables with MATLAB file format is not supported yet.');
                            case('TXT')
                                fileID = fopen(FluidProperties.RelPerm.TableFile{n,1}, 'r');
                                % Read lines from input file
                                Temp = textscan(fileID, '%s', 'Delimiter', '\n');
                                Temp = Temp{1};
                                Temp = Temp(2:end);
                                FluidProperties.RelPerm.TableData{n,1} = str2double(split(Temp,' '));
                                fclose(fileID);
                        end
                    end
            end
            
            % 8. Capillarity
            temp = strfind(obj.InputMatrix, 'CAPILLARITY');
            index = find(~cellfun('isempty', temp));
            FluidProperties.Capillarity.name = char(obj.InputMatrix(index + 1));
            if ~isempty(FluidProperties.Capillarity.name) 
                FluidProperties.Capillarity.wetting = str2double(obj.InputMatrix(index + 2));
            end
           
            % 9. Component properties
            temp = strfind(obj.InputMatrix, 'COMPONENT_PROPERTIES');
            index = find(~cellfun('isempty', temp));
            for i = 1:FluidProperties.NofComponents
                FluidProperties.ComponentProperties(i,:) = str2double(strsplit(char(obj.InputMatrix(index + i * 2))));
            end
            
            % 10. Gravity
            temp = strfind(obj.InputMatrix, 'GRAVITY');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
                FluidProperties.Gravity = 'OFF';
            else
                FluidProperties.Gravity = char(obj.InputMatrix(index + 1));
            end
        end
        function RockProperties = ReadRockProperties(obj)
            % 1. Density of the rock
            temp = strfind(obj.InputMatrix, 'ROCK_DENSITY');
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                RockProperties.Density = str2double(char(obj.InputMatrix{index+1}));
            else
                RockProperties.Density = 2750;
            end
            if isnan(RockProperties.Density)
                error('The rock density is not being read properly. Please check the input file!');
            end
            
            % 2. Compressiblity of the rock
            temp = strfind(obj.InputMatrix, 'ROCK_COMPRESSIBILITY');
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                RockProperties.Compressibility = str2double(char(obj.InputMatrix{index+1}));
            else
                RockProperties.Compressibility = 0;
            end
            if isnan(RockProperties.Compressibility)
                error('The rock compressibility is not being read properly. Please check the input file!');
            end
            
            % 3. Heat conductivity of the rock
            temp = strfind(obj.InputMatrix, 'ROCK_HEAT_CONDUCTIVITY');
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                RockProperties.HeatConductivity = str2double(char(obj.InputMatrix{index+1}));
            else
                RockProperties.HeatConductivity = 4;
            end
            if isnan(RockProperties.HeatConductivity)
                error('The rock heat conductivity is not being read properly. Please check the input file!');
            end
            
            % 4. Specific heat of the rock
            temp = strfind(obj.InputMatrix, 'ROCK_SPECIFIC_HEAT');
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                RockProperties.SpecificHeat = str2double(char(obj.InputMatrix{index+1}));
            else
                RockProperties.SpecificHeat = 790;
            end
            if isnan(RockProperties.SpecificHeat)
                error('The rock specific heat is not being read properly. Please check the input file!');
            end
        end
        function InitialConditions = ReadInitialConditions(obj)
            temp = strfind(obj.InputMatrix, 'ICLUDE_INITIAL_STATE');
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                
                InitialConditions.Include = 1;
                InitialConditions.File    = strcat(obj.Directory, '\', char(obj.InputMatrix(index+2)));
                fprintf('\n---> Reading "Initial State" file ...');
                InitialConditions.LoadedData      = load(InitialConditions.File);
                InitialConditions.LoadedData(:,1) = [];  % We do not need the first column as it is just the index of the cells
                InitialConditions.LoadedData(:,1) = InitialConditions.LoadedData(:,1) * 1e5;  % The stored data for pressure is in [Bar], we covert it to [Pa]
                fprintf(' ---> Completed.\n');
                
            else
                
                InitialConditions.Include = 0;
                InitialConditions.LoadedData = 0;
                
                temp = strfind(obj.InputMatrix, 'INITIAL_PRESSURE');
                index = find(~cellfun('isempty', temp));
                if ~isempty(index)
                    InitialConditions.Pressure = str2double(char(obj.InputMatrix(index + 1)));
                else
                    InitialConditions.Pressure = Nan;
                end
                
                temp = strfind(obj.InputMatrix, 'INITIAL_TEMPERATURE');
                index = find(~cellfun('isempty', temp));
                if ~isempty(index)
                    InitialConditions.Temperature = str2double(char(obj.InputMatrix(index + 1)));
                else
                    InitialConditions.Temperature = Nan;
                end
                
                temp = strfind(obj.InputMatrix, 'INITIAL_ENTHALPY');
                index = find(~cellfun('isempty', temp));
                if ~isempty(index)
                    InitialConditions.Enthalpy = str2double(char(obj.InputMatrix(index + 1)));
                else
                    InitialConditions.Enthalpy = Nan;
                end
                
                temp = strfind(obj.InputMatrix, 'INITIAL_SATURATION_1');
                index = find(~cellfun('isempty', temp));
                if ~isempty(index)
                    InitialConditions.Saturation_1 = str2double(char(obj.InputMatrix(index + 1)));
                else
                    InitialConditions.Saturation_1 = Nan;
                end
                
                temp = strfind(obj.InputMatrix, 'INITIAL_SATURATION_2');
                index = find(~cellfun('isempty', temp));
                if ~isempty(index)
                    InitialConditions.Saturation_2 = str2double(char(obj.InputMatrix(index + 1)));
                else
                    InitialConditions.Saturation_2 = Nan;
                end
                
            end
        end
        function WellsInfo = ReadWellsInfo(obj, SimulationInput)
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
                temp = regexp(WellInputMatrix, 'BOUNDARY_CONDITION', 'match');
                boundary_condition = find(~cellfun('isempty', temp));
                
                % Reading the coordinates of well trajectory
                Well.Coordinate.Type = WellInputMatrix{coordinate+1};
                switch Well.Coordinate.Type
                    case('IJK_LIST')
                        Well.Coordinate.Value = zeros(constraint-coordinate-2,3);
                        for p = 1 : size(Well.Coordinate.Value,1)
                            ijk = strsplit(WellInputMatrix{coordinate+2+p-1}, {'	',' ',',','[',']'} );
                            ijk = ijk(2:end-1);
                            ijk = strrep(ijk,'NX',num2str(SimulationInput.ReservoirProperties.Grid.N(1)));
                            ijk = strrep(ijk,'NY',num2str(SimulationInput.ReservoirProperties.Grid.N(2)));
                            ijk = strrep(ijk,'NZ',num2str(SimulationInput.ReservoirProperties.Grid.N(3)));
                            ijk = [str2double(ijk{1}), str2double(ijk{2}), str2double(ijk{3})];
                            Well.Coordinate.Value(p,:) = ijk;
                        end
                        if any(mod(Well.Coordinate.Value(:),1) ~= 0)
                            error('In well #%d, the ijk coordinates result in non-integer cell indeces. Check the input file!\n', w);
                        end
                    case('XYZ_LIST')
                        Well.Coordinate.Value = zeros(constraint-coordinate-2,3);
                        for p = 1 : size(Well.Coordinate.Value,1)
                            xyz = strsplit(WellInputMatrix{coordinate+2+p-1}, {'	',' ',',','[',']'} );
                            xyz = xyz(2:end-1);
                            xyz = strrep(xyz,'LX',num2str(SimulationInput.ReservoirProperties.size(1)));
                            xyz = strrep(xyz,'LY',num2str(SimulationInput.ReservoirProperties.size(2)));
                            xyz = strrep(xyz,'LZ',num2str(SimulationInput.ReservoirProperties.size(3)));
                            xyz = [str2double(xyz{1}), str2double(xyz{2}), str2double(xyz{3})];
                            Well.Coordinate.Value(p,:) = xyz;
                        end
                    case('CELL_INDEX_LIST')
                        Well.Coordinate.Value = [];
                        for p = 1 : constraint-coordinate-2
                            IndexList = strsplit(WellInputMatrix{coordinate+2+p-1}, {'	',' ',',',';','[',']'} );
                            IndexList = IndexList(2:end-1);
                            IndexList = strrep(IndexList,'NX',num2str(SimulationInput.ReservoirProperties.Grid.N(1)));
                            IndexList = strrep(IndexList,'NY',num2str(SimulationInput.ReservoirProperties.Grid.N(2)));
                            IndexList = strrep(IndexList,'NZ',num2str(SimulationInput.ReservoirProperties.Grid.N(3)));
                            IndexList = str2double(IndexList);
                            Well.Coordinate.Value = [ Well.Coordinate.Value, IndexList ];
                        end
                    otherwise
                        error('In well #%d, the coordination keyword "IJK_LIST" or "XYZ_LIST" or "CELL_INDEX_LIST" is missing. Check the input file!\n', w);
                end
                
                % Reading the contraint of well
                Well.Constraint.Name = char(WellInputMatrix(constraint+1));
                Well.Constraint.Value = str2double(WellInputMatrix(constraint+2));
                
                % Reading the formula type of well
                Well.Formula.Type = char(WellInputMatrix(formula+1));
                Well.Formula.Value = str2double(WellInputMatrix(formula+2));
                
                % Reading the boundary condition of the well (only for injection)
                switch SimulationInput.FluidProperties.FluidModel
                    case{'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                        Well.BoundaryCondition.Name  = char(       WellInputMatrix(boundary_condition+1) );
                        Well.BoundaryCondition.Value = str2double( WellInputMatrix(boundary_condition+2) );
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
            % 3.1 Number of reports --> report is written every
            % T_total/n_reports s; if 0 report at every dt
            index = find(strcmp(obj.SettingsMatrix, 'REPORTS'));
            SimulatorSettings.Reports = str2double(obj.SettingsMatrix(index+1));
            % 3.1 Number of previous reports so that the file indexing can continue from previous simulations
            index = find(strcmp(obj.SettingsMatrix, 'NUM_PREVIOUS_REPORTS'));
            if isempty(index)
                SimulatorSettings.NumOfPreviousReports = 0;
            else
                SimulatorSettings.NumOfPreviousReports = str2double(obj.SettingsMatrix(index+1));
            end
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
                    case('Geothermal_SinglePhase')
                        SimulatorSettings.Formulation = 'Geothermal_SinglePhase';
                    case('Geothermal_MultiPhase')
                        SimulatorSettings.Formulation = 'Geothermal_MultiPhase';
                    case('Molar')
                        SimulatorSettings.Formulation = 'Molar';
                    otherwise
                        error('The formulation defined in the input file is not valid');
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
                SimulatorSettings.ADMSettings.maxLevel = zeros(1+SimulationInput.FracturesProperties.NumOfFrac, 1);
                SimulatorSettings.ADMSettings.maxLevel(1) = str2double(obj.SettingsMatrix(x+1));
                SimulatorSettings.ADMSettings.CoarseningSwitch = zeros(1+SimulationInput.FracturesProperties.NumOfFrac, 1);
                SimulatorSettings.ADMSettings.CoarseningSwitch(1) = 1;
                SimulatorSettings.ADMSettings.Coarsening = zeros( 1+SimulationInput.FracturesProperties.NumOfFrac, 3, SimulatorSettings.ADMSettings.maxLevel(1) );
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
                SimulatorSettings.ADMSettings.Coupling = char(obj.SettingsMatrix(x+2));
                if ~strcmp(SimulatorSettings.ADMSettings.Coupling,'COUPLED') && ~strcmp(SimulatorSettings.ADMSettings.Coupling,'DECOUPLED')
                    SimulatorSettings.ADMSettings.Coupling = 'DECOUPLED';
                end
                for L = 1:SimulatorSettings.ADMSettings.maxLevel(1)
                    SimulatorSettings.ADMSettings.Coarsening(1,:,L) = [cx, cy, cz].^L; %Coarsening Factors: Cx1, Cy1; Cx2, Cy2; ...; Cxn, Cyn;
                end
                temp = strfind(obj.SettingsMatrix, 'PRESSURE_INTERPOLATOR');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.ADMSettings.PInterpolator = char(obj.SettingsMatrix(x+1));
                temp = strfind(obj.SettingsMatrix, 'HYPERBOLIC_INTERPOLATOR');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.ADMSettings.HInterpolator = char(obj.SettingsMatrix(x+1));

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
                    SimulatorSettings.ADMSettings.BF_MaxContrast = inf;
                else
                    SimulatorSettings.ADMSettings.BF_MaxContrast = str2double(obj.SettingsMatrix(Index+1));
                end
                
                % Maximum heterogeneity contrast of pEDFM alpha factors for basis functions
                temp = strfind(obj.SettingsMatrix, 'pEDFM_MAX_CONTRAST');
                Index = find(~cellfun('isempty', temp));
                if isempty(Index)
                    SimulatorSettings.ADMSettings.pEDFM_MaxContrast = inf;
                else
                    SimulatorSettings.ADMSettings.pEDFM_MaxContrast = str2double(obj.SettingsMatrix(Index+1));
                end
                
                % Correction of basis functions region of influence by cutting the values less than alpha  
                temp = strfind(obj.SettingsMatrix, 'BASIS_FUNCTION_REGION_OF_INFLUENCE_CUT');
                Index = find(~cellfun('isempty', temp));
                if isempty(Index)
                    SimulatorSettings.ADMSettings.BF_alpha = 0;
                else
                    SimulatorSettings.ADMSettings.BF_alpha = str2double(obj.SettingsMatrix(Index+1));
                end
                
                % ADM settings in the fractures
                if SimulationInput.FracturesProperties.isFractured
                    NumOfFrac = SimulationInput.FracturesProperties.NumOfFrac;
                    temp = strfind(obj.FractureMatrix, 'PROPERTIES');
                    frac_index = find(~cellfun('isempty', temp));
                    for f = 1 : NumOfFrac
                        frac_info_split = strsplit(obj.FractureMatrix{frac_index(f)},' ');
                        ADM_temp = regexprep(frac_info_split{9},' ' ,'');
                        ADM_temp = strsplit(ADM_temp, { '[' , ',' , ']' });
                        ADM_temp = [ str2double(ADM_temp(2)) , str2double(ADM_temp(3)) , str2double(ADM_temp(4)) , str2double(ADM_temp(5)) ];
                        SimulatorSettings.ADMSettings.CoarseningSwitch(1+f) = ADM_temp(1);
                        SimulatorSettings.ADMSettings.maxLevel(1+f) = ADM_temp(2);
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
                SimulatorSettings.MMsSettings.maxlevel = zeros(1+SimulationInput.FracturesProperties.NumOfFrac, 1);
                SimulatorSettings.MMsSettings.maxLevel(1) = str2double(obj.SettingsMatrix(x+1));
                SimulatorSettings.MMsSettings.CoarseningSwitch = zeros(1+SimulationInput.FracturesProperties.NumOfFrac, 1);
                SimulatorSettings.MMsSettings.CoarseningSwitch(1) = 1;
                SimulatorSettings.MMsSettings.Coarsening = zeros( 1+SimulationInput.FracturesProperties.NumOfFrac, 3, SimulatorSettings.MMsSettings.maxLevel(1) );
                temp = strfind(obj.SettingsMatrix, 'COARSENING_RATIOS');
                x = find(~cellfun('isempty', temp));
                cx = str2double(obj.SettingsMatrix(x+1));
                cy = str2double(obj.SettingsMatrix(x+2));
                cz = str2double(obj.SettingsMatrix(x+3));
                for L = 1:SimulatorSettings.MMsSettings.maxLevel(1)
                    SimulatorSettings.MMsSettings.Coarsening(1,:,L) = [cx, cy, cz].^L; %Coarsening Factors: Cx1, Cy1; Cx2, Cy2; ...; Cxn, Cyn;
                end
                
                temp = strfind(obj.SettingsMatrix, 'PRESSURE_INTERPOLATOR');
                x = find(~cellfun('isempty', temp));
                SimulatorSettings.MMsSettings.PInterpolator = char(obj.SettingsMatrix(x+1));
                
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
                    SimulatorSettings.MMsSettings.BF_MaxContrast = inf;
                else
                    SimulatorSettings.MMsSettings.BF_MaxContrast = str2double(obj.SettingsMatrix(Index+1));
                end
                
                % Maximum heterogeneity contrast of pEDFM alpha factors for basis functions
                temp = strfind(obj.SettingsMatrix, 'pEDFM_MAX_CONTRAST');
                Index = find(~cellfun('isempty', temp));
                if isempty(Index)
                    SimulatorSettings.MMsSettings.pEDFM_MaxContrast = inf;
                else
                    SimulatorSettings.MMsSettings.pEDFM_MaxContrast = str2double(obj.SettingsMatrix(Index+1));
                end
                
                % Correction of basis functions region of influence by cutting the values less than alpha  
                temp = strfind(obj.SettingsMatrix, 'BASIS_FUNCTION_REGION_OF_INFLUENCE_CUT');
                Index = find(~cellfun('isempty', temp));
                if isempty(Index)
                    SimulatorSettings.MMsSettings.BF_alpha = 0;
                else
                    SimulatorSettings.MMsSettings.BF_alpha = str2double(obj.SettingsMatrix(Index+1));
                end
                
                % MMs settings in fractures
                if SimulationInput.FracturesProperties.isFractured
                    NumOfFrac = SimulationInput.FracturesProperties.NumOfFrac;
                    temp = strfind(obj.FractureMatrix, 'PROPERTIES');
                    frac_index = find(~cellfun('isempty', temp));
                    for f = 1 : NumOfFrac
                        frac_info_split = strsplit(obj.FractureMatrix{frac_index(f)},' ');
                        MMs_temp = regexprep(frac_info_split{9},' ' ,'');
                        MMs_temp = strsplit(MMs_temp, { '[' , ',' , ']' });
                        MMs_temp = [ str2double(MMs_temp(2)) , str2double(MMs_temp(3)) , str2double(MMs_temp(4)) , str2double(MMs_temp(5)) ];
                        SimulatorSettings.MMsSettings.CoarseningSwitch(1+f) = MMs_temp(1);
                        SimulatorSettings.MMsSettings.maxLevel(1+f) = MMs_temp(2);
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
            
            SimulatorSettings.plotting.Software = char(obj.SettingsMatrix(xv+1)); %Matlab or ParaView/VTK
            SimulatorSettings.plotting.Format = char(obj.SettingsMatrix(xv+2)); % ASCII or BINARY
            
            temp = strfind(obj.SettingsMatrix, 'PLOT_INTERFACES');
            index = find(~cellfun('isempty', temp));
            if ~isempty(index)
                SimulatorSettings.plotting.PlotInterfaces = 1;
            else
                SimulatorSettings.plotting.PlotInterfaces = 0;
            end
            
            temp = strfind(obj.SettingsMatrix, 'PLOT_BASIS_FUNCTIONS');
            xv = find(~cellfun('isempty', temp));
            if ~isempty(xv)
                SimulatorSettings.PlotBasisFunctions = true;
            else
                SimulatorSettings.PlotBasisFunctions = false;
            end
        end
    end
end