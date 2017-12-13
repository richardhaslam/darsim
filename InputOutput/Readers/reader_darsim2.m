% Reader for Matteo's input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 December 2017
%Last modified: 13 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reader_darsim2 < reader
    properties
        %InputMatrix
        %SettingsMatrix
        %FractureMatrix
    end
    methods
        function obj = reader_darsim2(dir, file)
            obj@reader(dir, file);
        end
        function ReadInputFile(obj, Builder)
            %ReadInputFile
            fileID = fopen(obj.File, 'r');
            %// Read lines from input file
            matrix = textscan(fileID, '%s', 'Delimiter', '\n');
            obj.InputMatrix = matrix{1};
            fclose(fileID);
            %ReadSettingFile
            SettingsFile = strcat(obj.Directory, '/SimulatorSettings.txt');
            fileID = fopen(SettingsFile, 'r');
            matrix = textscan(fileID, '%s', 'Delimiter', '\n');
            obj.SettingsMatrix = matrix{1};
            fclose(fileID);
            %ReadSettingFile
            temp = strfind(obj.InputMatrix, 'FRACTURED');
            fractured = find(~cellfun('isempty', temp));
            if isempty(fractured)
                fractured = 0;
            else
                fractured = 1;
            end
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
            SimulationInput.TotalTime = str2double(obj.InputMatrix(xv + 1))*24*3600;
            
            % 3. Read reservoir properties
            SimulationInput.ReservoirProperties = obj.ReadReservoirProperties();
            
            % 4. Read fluid properties
            SimulationInput.FluidProperties = obj.ReadFluidProperties();
                                  
            % 5. Initial conditions
            temp = strfind(obj.InputMatrix, 'INIT');
            index = find(~cellfun('isempty', temp));
            SimulationInput.Init = str2double(strsplit(char(obj.InputMatrix(index + 1))));
            
            % 6. Read wells info
            SimulationInput.WellsInfo = obj.ReadWellsInfo(SimulationInput);
            
            % 7. Read Fractures' Properties            
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
            % 1. size of the reservoir
            temp = strfind(obj.InputMatrix, 'DIMENS'); % Search a specific string and find all rows containing matches
            index = find(~cellfun('isempty', temp));
            ReservoirProperties.size = [str2double(obj.InputMatrix{index+1});...
                                        str2double(obj.InputMatrix{index+2});
                                        str2double(obj.InputMatrix{index+3})];
            % 2. GridSize
            temp = strfind(obj.InputMatrix, 'SPECGRID');
            index = find(~cellfun('isempty', temp));
            ReservoirProperties.Grid.N = [str2double(obj.InputMatrix{index+1});...
                                          str2double(obj.InputMatrix{index+2});
                                          str2double(obj.InputMatrix{index+3})];
            
            % 3. Permeability
            perm = zeros(3, 1);
            temp = strfind(obj.InputMatrix, 'PERMX');
            perm(1) = find(~cellfun('isempty', temp));
            temp = strfind(obj.InputMatrix, 'PERMY');
            perm(2) = find(~cellfun('isempty', temp));
            temp = strfind(obj.InputMatrix, 'PERMZ');
            perm(3) = find(~cellfun('isempty', temp));
            for i=1:3
                if strcmp(obj.InputMatrix(perm(i) - 1), 'INCLUDE')
                    ReservoirProperties.PermInclude(i) = 1;
                    ReservoirProperties.PermFile(i) = strcat('../Permeability/', char(obj.InputMatrix(perm(i) +1)));
                else
                    ReservoirProperties.PermInclude(i) = 0;
                    ReservoirProperties.Perm(i) = str2double(obj.InputMatrix(perm(i) +1));
                end
            end

            % 4. Porosity 
            temp = strfind(obj.InputMatrix, 'POR');
            index = find(~cellfun('isempty', temp));
            ReservoirProperties.phi = str2double(obj.InputMatrix(index + 1));
            
            % 5. Temperature
            temp = strfind(obj.InputMatrix, 'TEMPERATURE (K)');
            index = find(~cellfun('isempty', temp));
            ReservoirProperties.Temperature = str2double(obj.InputMatrix(index + 1));;
        end
        function FluidProperties = ReadFluidProperties(obj)
            %%%%%%%%%%%%%FLUID PROPERTIES%%%%%%%%%%%%%%%%
            % 1. Fluid model
            temp = strfind(obj.InputMatrix, 'FLUID MODEL');
            index = find(~cellfun('isempty', temp));
            FluidProperties.FluidModel = str2double(obj.InputMatrix{index+1});
            FluidProperties.NofPhases = str2double(obj.InputMatrix{index+3});
            FluidProperties.NofComponents = str2double(obj.InputMatrix{index+5});
            % 2. Density
            temp = strfind(obj.InputMatrix, 'DENSITY');
            index_density = find(~cellfun('isempty', temp));
            % 3. Viscosity
            temp = strfind(obj.InputMatrix, 'VISCOSITY');
            index_viscosity = find(~cellfun('isempty', temp));
            % 4. Compressibility
            temp = strfind(obj.InputMatrix, 'COMPRESSIBILITY');
            index_comp = find(~cellfun('isempty', temp));
            for i=1:FluidProperties.NofPhases
                FluidProperties.Density(i) = str2double(obj.InputMatrix{index_density+(i*2)});
                FluidProperties.mu(i) = str2double(obj.InputMatrix{index_viscosity+(i*2)});
                FluidProperties.Comp(i) = str2double(obj.InputMatrix{index_comp+(i*2)});
            end
            % 5. Relative permeability
            temp = strfind(obj.InputMatrix, 'RELPERM');
            index = find(~cellfun('isempty', temp));
            FluidProperties.RelPerm.name = char(obj.InputMatrix(index + 1));
            FluidProperties.RelPerm.s_irr(1) = str2double(obj.InputMatrix(index + 3));
            FluidProperties.RelPerm.s_irr(2) = str2double(obj.InputMatrix(index + 5));
            FluidProperties.RelPerm.Kre(1) = str2double(obj.InputMatrix(index + 7));
            FluidProperties.RelPerm.Kre(2) = str2double(obj.InputMatrix(index + 9));
            FluidProperties.RelPerm.n(1) = str2double(obj.InputMatrix(index + 11));
            FluidProperties.RelPerm.n(2) = str2double(obj.InputMatrix(index + 13));

            % 6. Capillarity
            temp = strfind(obj.InputMatrix, 'CAPILLARITY');
            index = find(~cellfun('isempty', temp));
            FluidProperties.Capillarity.name = char(obj.InputMatrix(index + 1));
            if ~isempty(FluidProperties.Capillarity.name) 
                FluidProperties.Capillarity.wetting = str2double(obj.InputMatrix(index + 2));
            end
           
            % 7. Component properties
            temp = strfind(obj.InputMatrix, 'COMPONENT PROPERTIES');
            index = find(~cellfun('isempty', temp));
            %FluidProperties.ComponentProperties;
            
            % 9. Gravity
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
            temp = regexp(obj.InputMatrix, 'INJ\d', 'match');
            inj = find(~cellfun('isempty', temp));
            WellsInfo.NofInj = length(inj);
            for i=1:WellsInfo.NofInj 
                WellsInfo.Inj(i).Coord = obj.ReadCoordinates(inj(i), SimulationInput);
                WellsInfo.Inj(i).Constraint.name = char(obj.InputMatrix(inj(i) + 7));
                WellsInfo.Inj(i).Constraint.value = str2double(obj.InputMatrix(inj(i) + 8));
                WellsInfo.Inj(i).PI.type = char(obj.InputMatrix(inj(i) + 9));
                WellsInfo.Inj(i).PI.value = str2double(obj.InputMatrix(inj(i) + 10));
            end
            
            temp = regexp(obj.InputMatrix, 'PROD\d', 'match');
            prod = find(~cellfun('isempty', temp));
            WellsInfo.NofProd = length(prod);
            for i=1:WellsInfo.NofProd
                WellsInfo.Prod(i).Coord = obj.ReadCoordinates(prod(i), SimulationInput);
                WellsInfo.Prod(i).Constraint.name = char(obj.InputMatrix(prod(i) + 7));
                WellsInfo.Prod(i).Constraint.value = str2double(obj.InputMatrix(prod(i) + 8));
                WellsInfo.Prod(i).PI.type = char(obj.InputMatrix(prod(i) + 9));
                WellsInfo.Prod(i).PI.value = str2double(obj.InputMatrix(prod(i) + 10));
            end
        end
        function Coord = ReadCoordinates(obj, index, SimulationInput)
            % Read coordinates of the wells
            %i_init
            Well_Coord_Temp = strsplit(obj.InputMatrix{index + 1}, ' ');
            if sum( strcmp('NX' , Well_Coord_Temp) ) > 0
                Well_Coord_Temp = strsplit(obj.InputMatrix{index + 1}, ' ');
                if length(Well_Coord_Temp)>1
                    if Well_Coord_Temp{2}=='-',  i_init =SimulationInput.ReservoirProperties.Grid.N(1)-1;
                    else,  error('For Injection Well Coordination, you can only use "NX" with minus sign "-"!');
                    end
                else
                    i_init =SimulationInput.ReservoirProperties.Grid.N(1);
                end
            else
                i_init = str2double(obj.InputMatrix{index + 1});
            end
            %i_final
            Well_Coord_Temp = strsplit(obj.InputMatrix{index + 2}, ' ');
            if sum( strcmp('NX' , Well_Coord_Temp) ) > 0
                Well_Coord_Temp = strsplit(obj.InputMatrix{index + 2}, ' ');
                if length(Well_Coord_Temp)>1
                    if Well_Coord_Temp{2}=='-',  i_final =SimulationInput.ReservoirProperties.Grid.N(1)-1;
                    else,  error('For Injection Well Coordination, you can only use "NX" with minus sign "-"!');
                    end
                else
                    i_final =SimulationInput.ReservoirProperties.Grid.N(1);
                end
            else
                i_final = str2double(obj.InputMatrix{index + 2});
            end
            %j_init
            Well_Coord_Temp = strsplit(obj.InputMatrix{index + 3}, ' ');
            if sum( strcmp('NY' , Well_Coord_Temp) ) > 0
                Well_Coord_Temp = strsplit(obj.InputMatrix{index + 3}, ' ');
                if length(Well_Coord_Temp)>1
                    if Well_Coord_Temp{2}=='-',  j_init =SimulationInput.ReservoirProperties.Grid.N(2)-1;
                    else,  error('For Injection Well Coordination, you can only use "NY" with minus sign "-"!');
                    end
                else
                    j_init =SimulationInput.ReservoirProperties.Grid.N(2);
                end
            else
                j_init = str2double(obj.InputMatrix{index + 3});
            end
            %j_final
            Well_Coord_Temp = strsplit(obj.InputMatrix{index + 4}, ' ');
            if sum( strcmp('NY' , Well_Coord_Temp) ) > 0
                Well_Coord_Temp = strsplit(obj.InputMatrix{index + 4}, ' ');
                if length(Well_Coord_Temp)>1
                    if Well_Coord_Temp{2}=='-',  j_final =SimulationInput.ReservoirProperties.Grid.N(2)-1;
                    else,  error('For Injection Well Coordination, while you can only use "NY" with minus sign "-"!');
                    end
                else
                    j_final =SimulationInput.ReservoirProperties.Grid.N(2);
                end
            else
                j_final = str2double(obj.InputMatrix{index + 4});
            end
            %k_init
            Well_Coord_Temp = strsplit(obj.InputMatrix{index + 5}, ' ');
            if sum( strcmp('NZ' , Well_Coord_Temp) ) > 0
                Well_Coord_Temp = strsplit(obj.InputMatrix{index + 5}, ' ');
                if length(Well_Coord_Temp)>1
                    if Well_Coord_Temp{2}=='-',  k_init =SimulationInput.ReservoirProperties.Grid.N(3)-1;
                    else,  error('For Injection Well Coordination, you can only use "NZ" with minus sign "-"!');
                    end
                else
                    k_init =SimulationInput.ReservoirProperties.Grid.N(3);
                end
            else
                k_init = str2double(obj.InputMatrix{index + 5});
            end
            Well_Coord_Temp = strsplit(obj.InputMatrix{index + 6}, ' ');
            %k_final
            if sum( strcmp('NZ' , Well_Coord_Temp) ) > 0
                Well_Coord_Temp = strsplit(obj.InputMatrix{index + 6}, ' ');
                if length(Well_Coord_Temp)>1
                    if Well_Coord_Temp{2}=='-',  k_final =SimulationInput.ReservoirProperties.Grid.N(3)-1;
                    else,  error('For Injection Well Coordination, while you can only use "NZ" with minus sign "-"!');
                    end
                else
                    k_final =SimulationInput.ReservoirProperties.Grid.N(3);
                end
            else
                k_final = str2double(obj.InputMatrix{index + 6});
            end
            Coord = [i_init, i_final; j_init, j_final; k_init, k_final];
        end
        function SimulatorSettings = ReadSimulatorSettings(obj, SimulationInput)
            % 1. Max number of time-steps
            temp = strfind(obj.SettingsMatrix, 'TIMESTEPS');
            xv = find(~cellfun('isempty', temp));
            SimulatorSettings.MaxNumTimeSteps = str2double(obj.SettingsMatrix(xv+1));
            % 2. Min and max time-step sizes
            temp = strfind(obj.SettingsMatrix, 'MINMAXDT');
            xv = find(~cellfun('isempty', temp));
            if ~isempty(xv)
                SimulatorSettings.MinMaxdt = str2double(strsplit(char(obj.SettingsMatrix(xv+1))));
            else
                % default value
                SimulatorSettings.MinMaxdt = [0.1, 30];
            end
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
                SimulatorSettings.Tolerance = str2double(obj.SettingsMatrix(index + 2));
            else
                SimulatorSettings.CouplingType ='Sequential';
                temp = strfind(obj.SettingsMatrix, 'SEQUENTIAL'); 
                index = find(~cellfun('isempty', temp));
                SimulatorSettings.MaxIterations = str2double(obj.SettingsMatrix(index + 1)); 
                SimulatorSettings.Tolerance     = str2double(obj.SettingsMatrix(index + 2));
                SimulatorSettings.cfl           = str2double(obj.SettingsMatrix(index + 3)); 
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
            
            %% ADM settings
            temp = strfind(obj.SettingsMatrix, 'ADM');
            adm = find(~cellfun('isempty', temp));
            if str2double(obj.SettingsMatrix(adm + 1)) == 1
                SimulatorSettings.ADM = 'active';
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
                if isempty(SimulatorSettings.ADMSettings.maxLevel(1)) || isempty(SimulatorSettings.ADMSettings.Coarsening(1,:,:)) || isempty(SimulatorSettings.ADMSettings.key) || isempty(SimulatorSettings.ADMSettings.tol) || isempty(SimulatorSettings.ADMSettings.PInterpolator)
                    error('DARSIM2 ERROR: Missing ADM settings! Povide LEVELS, COARSENING_CRITERION, COARSENING_RATIOS, TOLERANCE, PRESSURE_INTERPOLATOR');
                end
                
                temp = strfind(obj.SettingsMatrix, 'COUPLED');
                temp = find(~cellfun('isempty', temp));
                if isempty(temp)
                    SimulatorSettings.ADMSettings.BFtype = 'decoupled';
                else
                    SimulatorSettings.ADMSettings.BFtype = 'coupled';
                end
                
                % ADM settings in the fractures
                if SimulationInput.FracturesProperties.Fractured
                    for f = 1 : SimulationInput.FracturesProperties.NrOfFrac
                        frac_info_split = strsplit(obj.FractureMatrix{1}{frac_index(f)},' ');
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
                                SimulatorSettings.ADMSettings.Coarsening(1+f,:,L) = [ADM_temp(3), ADM_temp(4), 1].^obj.ADMSettin.maxLevel(1+f);
                            end
                        end
                    end
                end
            end
            
            %% Multilevel Multiscale settings
            temp = strfind(obj.SettingsMatrix, 'MMs');
            mms = find(~cellfun('isempty', temp));
            if str2double(obj.SettingsMatrix(mms + 1)) == 1
                SimulatorSettings.MMs = 'active';
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
                
                temp = strfind(obj.SettingsMatrix, 'MSFE');
                temp = find(~cellfun('isempty', temp));
                if isempty(temp)
                    SimulatorSettings.MMsSettings.MSFE = 0;
                else
                    SimulatorSettings.MMsSettings.MSFE = 1;
                end
                
                temp = strfind(obj.SettingsMatrix, 'COUPLED');
                temp = find(~cellfun('isempty', temp));
                if isempty(temp)
                    SimulatorSettings.MMsSettings.BFtype = 'DECOUPLED';
                else
                    SimulatorSettings.MMsSettings.BFtype = 'COUPLED';
                end
                
                % MMs settings in fractures
                if obj.Fractured
                    for f = 1 : SimulationInput.FracturesProperties.NrOfFrac
                        frac_info_split = strsplit(obj.FractureMatrix{1}{frac_index(f)},' ');
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
            
            %%%%%%%%%%%%%OPTIONS%%%%%%%%%%%%%%%%
            temp = strfind(obj.SettingsMatrix, 'OUTPUT'); 
            xv = find(~cellfun('isempty', temp));
            SimulatorSettings.plotting = char(obj.SettingsMatrix(xv+1)); %Matlab or VTK
        end
    end
end