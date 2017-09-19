 % Builder Builds all objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 24 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef builder < handle
    properties
        ProblemName
        TotalTime
        size
        temperature
        grid
        Gravity
        perm
        pert
        por
        density
        viscosity
        relperm
        compressibility
        capillarity
        foam
        Comp_Type
        Comp_Prop
        Init
        inj
        prod
        MaxNumTimeSteps
        reports
        CouplingType
        coupling
        transport
        plotting
        adm
        flash
        ADM
        LinearSolver
        Formulation = 'Natural';
        StopCriterion = 'MAX TIME';
        Fractured = 0;
        incompressible
    end
    methods
        function FindKeyWords(obj, inputMatrix, SettingsMatrix)
            %%%%%% Name of the Problem %%%%%%
            temp = strfind(inputMatrix{1}, 'TITLE'); % Search a specific string and find all rows containing matches
            obj.ProblemName = char(inputMatrix{1}(find(~cellfun('isempty', temp)) + 1));
            
            %%%%%%%%%%%%%PROPERTIES OF THE RESERVOIR%%%%%%%%%%%%%%%%
            temp = strfind(inputMatrix{1}, 'DIMENS'); % Search a specific string and find all rows containing matches
            obj.size = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'SPECGRID');
            obj.grid = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'PERMX');
            obj.perm = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'PERTURB');
            obj.pert = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'POR');
            obj.por = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'TEMPERATURE (K)');
            obj.temperature = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'TOTALTIME');
            xv = find(~cellfun('isempty', temp));
            obj.TotalTime = str2double(inputMatrix{1}(xv + 1))*24*3600;
            
            %%%%%%%%%%%%%FLUID PROPERTIES%%%%%%%%%%%%%%%%
            temp = strfind(inputMatrix{1}, 'DENSITY'); 
            obj.density = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'VISCOSITY');
            obj.viscosity = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'RELPERM');
            obj.relperm = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'COMPRESSIBILITY');
            obj.compressibility = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'CAPILLARITY');
            obj.capillarity = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'FOAM');
            obj.foam = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'FLUID MODEL');
            obj.Comp_Type = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'COMPONENT PROPERTIES');
            obj.Comp_Prop = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'GRAVITY');
            index = find(~cellfun('isempty', temp));
            if isempty(index)
               obj.Gravity = 'OFF'; 
            else
                obj.Gravity = char(inputMatrix{1}(index + 1));
            end
            
            %%%%%%%%%%%%%INITIAL CONDITIONS%%%%%%%%%%%%%%%%
            temp = strfind(inputMatrix{1}, 'INIT');
            init = find(~cellfun('isempty', temp));
            obj.Init = str2double(strsplit(char(inputMatrix{1}(init + 1))));
            
            %%%%%%%%%%%%%WELLS%%%%%%%%%%%%%%%%
            temp = regexp(inputMatrix{1}, 'INJ\d', 'match');
            obj.inj = find(~cellfun('isempty', temp));
            temp = regexp(inputMatrix{1}, 'PROD\d', 'match');
            obj.prod = find(~cellfun('isempty', temp));
            
            %%%%%%%%%%%%%%Fractures' Properties%%%%%%%%%%%%%%%%
            % Check the main input file and look for FRACTURED
            temp = strfind(inputMatrix{1}, 'FRACTURED');
            obj.Fractured = find(~cellfun('isempty', temp));
            if isempty(obj.Fractured)
                obj.Fractured = 0;
            else
                obj.Fractured = 1;
            end
           
            %%%%%%%%%%%%%%%SIMULATOR'S SETTINGS%%%%%%%%%%%
            temp = strfind(SettingsMatrix{1}, 'TIMESTEPS');
            xv = find(~cellfun('isempty', temp));
            obj.MaxNumTimeSteps = str2double(SettingsMatrix{1}(xv+1));
            temp = strfind(SettingsMatrix{1}, 'REPORTS');
            xv = find(~cellfun('isempty', temp));
            obj.reports = str2double(SettingsMatrix{1}(xv+1));
            temp = strfind(SettingsMatrix{1}, 'FIM'); 
            obj.coupling = find(~cellfun('isempty', temp));
            temp = strfind(SettingsMatrix{1}, 'FORMULATION'); 
            xv = find(~cellfun('isempty', temp));
            obj.Formulation =  char(SettingsMatrix{1}(xv+1));
            temp = strfind(SettingsMatrix{1}, 'LINEARSOLVER'); 
            xv = find(~cellfun('isempty', temp));
            obj.LinearSolver =  char(SettingsMatrix{1}(xv+1));
            if obj.coupling ~= 0
                obj.CouplingType = 'FIM';
            else
                obj.CouplingType ='Sequential';
                temp = strfind(SettingsMatrix{1}, 'SEQUENTIAL'); 
                obj.coupling = find(~cellfun('isempty', temp));
                temp = strfind(SettingsMatrix{1}, 'IMPSAT');
                obj.transport = find(~cellfun('isempty', temp));
                obj.Formulation = 'Sequential';
            end
            
            temp = strfind(SettingsMatrix{1}, 'ADM');
            obj.adm = find(~cellfun('isempty', temp));
            
            %%%%%%%%%%%%%OPTIONS%%%%%%%%%%%%%%%%
            temp = strfind(SettingsMatrix{1}, 'OUTPUT'); 
            xv = find(~cellfun('isempty', temp));
            obj.plotting = char(SettingsMatrix{1}(xv+1)); %Matlab or VTK
            
        end
        function simulation = BuildSimulation(obj, inputMatrix, SettingsMatrix, FractureMatrix)
            simulation = Reservoir_Simulation();
            simulation.DiscretizationModel = obj.BuildDiscretization(inputMatrix, FractureMatrix, SettingsMatrix);
            simulation.ProductionSystem = obj.BuildProductionSystem(inputMatrix, FractureMatrix, simulation.DiscretizationModel);            
            simulation.FluidModel = obj.BuildFluidModel(inputMatrix, simulation.ProductionSystem);
            simulation.Formulation = obj.BuildFormulation(simulation.DiscretizationModel, simulation.FluidModel, simulation.ProductionSystem);
            simulation.TimeDriver = obj.BuildTimeDriver(SettingsMatrix);
            simulation.Summary = obj.BuildSummary(simulation);
            
            % Define Properties
            obj.DefineProperties(simulation.ProductionSystem, simulation.FluidModel, simulation.DiscretizationModel);
            
            %% Define Initialization procedure
            N = simulation.DiscretizationModel.ReservoirGrid.N;
            VarValues = ones(N, length(obj.Init));
            for i=1:length(obj.Init)
                VarValues(:, i) = VarValues(:, i) * obj.Init(i);
            end
            switch(simulation.FluidModel.name)
                case('SinglePhase') 
                    VarNames = {'P_1', 'S_1'};
                    VarValues(:, 2) = 1;
                    simulation.Initializer = initializer_singlephase(VarNames, VarValues);
                case('Immiscible')
                    VarNames = {'P_2', 'S_1', 'S_2'};
                    %index = 1;
                    %VarValues(index, 2) = 1;
                    simulation.Initializer = initializer(VarNames, VarValues);
                otherwise
                    VarNames = {'P_2', 'z_1', 'z_2'};
                    simulation.Initializer = initializer_hydrostatic(VarNames, VarValues);
            end
        end
        function Discretization = BuildDiscretization(obj, inputMatrix, FractureMatrix, SettingsMatrix)
            %% 1. Create fine-scale grids
            % 1a. Reservoir Grid
            nx = str2double(inputMatrix(obj.grid + 1));
            ny = str2double(inputMatrix(obj.grid + 2));
            nz = str2double(inputMatrix(obj.grid + 3));
            Nm = nx*ny*nz;
            ReservoirGrid = cartesian_grid(nx, ny, nz);
            NrOfFrac = 0;
            % 1b. Fractures Grid
            if obj.Fractured
                temp = strfind(FractureMatrix{1}, 'NUM_FRACS');
                index = find(~cellfun('isempty', temp));
                temp = strsplit(FractureMatrix{1}{index},' ');
                NrOfFrac = str2double( temp{end} );
                FracturesGrid = fractures_grid(NrOfFrac);
                temp = strfind(FractureMatrix{1}, 'PROPERTIES');
                frac_index = find(~cellfun('isempty', temp));
                
                Nf = zeros(NrOfFrac,1);
                for f = 1 : NrOfFrac
                    % Creat cartesian grid in each fracture
                    frac_info_split = strsplit(FractureMatrix{1}{frac_index(f)},' ');
                    grid_temp = strsplit(frac_info_split{8}, 'x');
                    nx = str2double(grid_temp{1});
                    ny = str2double(grid_temp{2});
                    nz = 1;
                    FractureGrid = cartesian_grid(nx, ny, nz);
                    
                    % Add fracture grid coordinates (it's for plotting purposes) 
                    temp = strfind(FractureMatrix{1}, 'GRID_COORDS_X');
                    frac_grid_coords_x = find(~cellfun('isempty', temp));
                    temp = strfind(FractureMatrix{1}, 'GRID_COORDS_Y');
                    frac_grid_coords_y = find(~cellfun('isempty', temp));
                    temp = strfind(FractureMatrix{1}, 'GRID_COORDS_Z');
                    frac_grid_coords_z = find(~cellfun('isempty', temp));
                    frac_grid_coords_x_split = strsplit(FractureMatrix{1}{frac_grid_coords_x(f)},' ');
                    frac_grid_coords_y_split = strsplit(FractureMatrix{1}{frac_grid_coords_y(f)},' ');
                    frac_grid_coords_z_split = strsplit(FractureMatrix{1}{frac_grid_coords_z(f)},' ');
                    frac_grid_coords_x_split = str2double(frac_grid_coords_x_split);
                    frac_grid_coords_y_split = str2double(frac_grid_coords_y_split);
                    frac_grid_coords_z_split = str2double(frac_grid_coords_z_split);
                    frac_grid_coords_x_split(1) = [];
                    frac_grid_coords_y_split(1) = [];
                    frac_grid_coords_z_split(1) = [];
                    FractureGrid.GridCoords = [ frac_grid_coords_x_split' , frac_grid_coords_y_split' , frac_grid_coords_z_split' ];
                    FracturesGrid.AddGrid(FractureGrid, f);
                    Nf(f) = nx*ny*nz;
                end
                
                % Reading the non-neighboring connectivities of frac-frac and frac-matrix
                CrossConnections = cross_connections();
                temp = strfind(FractureMatrix{1}, 'FRACCELL');
                frac_cell_index = find(~cellfun('isempty', temp));
                
                temp = strfind(FractureMatrix{1}, 'ROCK_CONN');
                frac_rockConn_index = find(~cellfun('isempty', temp));
                
                temp = strfind(FractureMatrix{1}, 'FRAC_CONN');
                frac_fracConn_index = find(~cellfun('isempty', temp));
                
                n_phases = str2double(inputMatrix(obj.Comp_Type + 3)); % Number of phases (useful to define size of some objects)
                for f = 1 : NrOfFrac
                    % looping over all global fracture cells
                    for If = 1:length(frac_cell_index)
                        fracCell_info_split = strsplit(FractureMatrix{1}{frac_cell_index(If)},{' ','	'});
                        
                        % frac-marix conn
                        temp = frac_rockConn_index - frac_cell_index(If); temp(temp<0) = max(temp) +1;
                        [~ , frac_rockConn_index_start] = min(temp);
                        for Im = 1:str2double(fracCell_info_split{3})
                            frac_rockConn_info_split = strsplit(FractureMatrix{1}{frac_rockConn_index(frac_rockConn_index_start+Im-1)},{' ','	'});
                            CrossConnections(If,1).Cells(Im,1) = str2double(frac_rockConn_info_split{2})+1;
                            CrossConnections(If,1).T_Geo(Im,1) = str2double(frac_rockConn_info_split{3});
                        end
                        
                        % frac-frac conn
                        temp = frac_fracConn_index - frac_cell_index(If); temp(temp<0) = max(temp)+1;
                        [~ , frac_fracConn_index_start] = min(temp);
                        
                        Counter = 1;
                        for Ig = 1:str2double(fracCell_info_split{5})
                            frac_fracConn_info_split = strsplit(FractureMatrix{1}{frac_fracConn_index(frac_fracConn_index_start+Ig-1)},{' ','	'});
                            If_Other_Global = Nm + sum( Nf( 1 : str2double(frac_fracConn_info_split{2}) +1-1 ) ) + str2double(frac_fracConn_info_split{3})+1;
                            if If_Other_Global > If + Nm
                                CrossConnections(If,1).Cells(Im+Counter,1) = If_Other_Global;
                                CrossConnections(If,1).T_Geo(Im+Counter,1) = str2double(frac_fracConn_info_split{4});
                                Counter = Counter + 1;
                            end
                        end
                        CrossConnections(If,1).UpWind = zeros(length(CrossConnections(If,1).Cells), n_phases);
                        CrossConnections(If,1).U_Geo = zeros(length(CrossConnections(If,1).Cells), n_phases);
                    end
                end
            end
            
            %% 2. Define your discretization Model (choose between FS and ADM)
            if (str2double(SettingsMatrix(obj.adm + 1)) == 0 )
                % Fine-scale discretization model
                obj.ADM = 'inactive';
                Discretization = FS_Discretization_model();
            else
                % ADM discretization model
                obj.ADM = 'active';
                maxLevel = ones(1+NrOfFrac , 1);
                % ADM grid for reservoir
                temp = strfind(SettingsMatrix, 'LEVELS');
                x = find(~cellfun('isempty', temp));
                maxLevel(1) = str2double(SettingsMatrix(x+1));
                Coarsening = zeros( 1+NrOfFrac, 3, maxLevel(1) );
                temp = strfind(SettingsMatrix, 'COARSENING_RATIOS');
                x = find(~cellfun('isempty', temp));
                cx = str2double(SettingsMatrix(x+1));
                cy = str2double(SettingsMatrix(x+2));
                cz = str2double(SettingsMatrix(x+3));
                temp = strfind(SettingsMatrix, 'COARSENING_CRITERION');
                x = find(~cellfun('isempty', temp));
                gridselcriterion = char(SettingsMatrix(x+1));
                temp = strfind(SettingsMatrix, 'VARIABLE');
                x = find(~cellfun('isempty', temp));
                key = char(SettingsMatrix(x+1));
                temp = strfind(SettingsMatrix, 'TOLERANCE');
                x = find(~cellfun('isempty', temp));
                tol = str2double(SettingsMatrix(x+1));
                for L = 1:maxLevel(1)
                    Coarsening(1,:,L) = [cx, cy, cz].^L; %Coarsening Factors: Cx1, Cy1; Cx2, Cy2; ...; Cxn, Cyn;
                end
                temp = strfind(SettingsMatrix, 'PRESSURE_INTERPOLATOR');
                x = find(~cellfun('isempty', temp));
                if isempty(maxLevel(1)) || isempty(Coarsening(1,:,:)) || isempty(key) || isempty(tol) || isempty(x)
                    error('DARSIM2 ERROR: Missing ADM settings! Povide LEVELS, COARSENING_CRITERION, COARSENING_RATIOS, TOLERANCE, PRESSURE_INTERPOLATOR');
                end 
                % ADM grid for fractures
                if obj.Fractured
                    for f = 1 : NrOfFrac
                        frac_info_split = strsplit(FractureMatrix{1}{frac_index(f)},' ');
                        ADM_temp = regexprep(frac_info_split{9},' ' ,'');
                        ADM_temp = strsplit(ADM_temp, { '[' , ',' , ']' });
                        ADM_temp = [ str2double(ADM_temp(2)) , str2double(ADM_temp(3)) , str2double(ADM_temp(4)) , str2double(ADM_temp(5)) ];
                        if ADM_temp(1)
                            maxLevel(1+f) = ADM_temp(2);
                        else
                            maxLevel(1+f) = 0;
                        end
                        for L = 1:maxLevel(1)
                            if L <= maxLevel(1+f)
                                Coarsening(1+f,:,L) = [ADM_temp(3), ADM_temp(4), 1].^L;
                            else
                                Coarsening(1+f,:,L) = [ADM_temp(3), ADM_temp(4), 1].^maxLevel(1+f);
                            end
                        end
                    end
                end
                % Create the operatorshandler
                operatorshandler = operators_handler(Coarsening(1,:,:));
                % a.1 Pressure prolongation builder
                switch (char(SettingsMatrix(x+1)))
                    case ('Constant')
                        prolongationbuilder = prolongation_builder_constant(maxLevel(1));
                    otherwise
                        prolongationbuilder = prolongation_builder_MSPressure(maxLevel(1), Coarsening(:,:,1) );
                        if ~obj.Fractured
                            prolongationbuilder.BFUpdater = bf_updater_ms();
                        else
                            prolongationbuilder.BFUpdater = bf_updater_FAMS();
                        end
                        if strcmp(char(SettingsMatrix(x+1)), 'Homogeneous')
                            prolongationbuilder.BFUpdater.MaxContrast = 1;
                        else
                            prolongationbuilder.BFUpdater.MaxContrast = 10^-2;
                        end
                end
                operatorshandler.AddProlongationBuilder(prolongationbuilder, 1);
                % a.2 Hyperbolic variables operators builder
                n_phases = str2double(inputMatrix(obj.Comp_Type + 3));
                test = 'Constant';
                for i = 2:n_phases
                    switch(test)
                        case('Constant')
                            prolongationbuilder = prolongation_builder_constant(maxLevel(1));
                        case('MultiscaleSat')
                            prolongationbuilder = prolongation_builder_MSHyperbolic(maxLevel(1));
                    end
                    operatorshandler.AddProlongationBuilder(prolongationbuilder, i);
                end
                
                % b. Grid selection criterion (time\space based)
                switch (gridselcriterion)
                    case('dfdx')
                        gridselector = adm_grid_selector_delta(tol, key);
                    case('dfdt')
                        gridselector = adm_grid_selector_time(tol, key);
                end
                Discretization = ADM_Discretization_model(maxLevel, Coarsening);
                Discretization.AddADMGridSelector(gridselector);
                Discretization.AddOperatorsHandler(operatorshandler);
            end
            %% 3. Add Grids to the Discretization Model
            Discretization.AddReservoirGrid(ReservoirGrid);
            if obj.Fractured
                Discretization.AddFracturesGrid(FracturesGrid);
                Discretization.AddCrossConnections(CrossConnections);
            end
        end
        function ProductionSystem = BuildProductionSystem (obj, inputMatrix, FractureMatrix, DiscretizationModel)
            ProductionSystem = Production_System();
            %% RESERVOIR
            Lx = str2double(inputMatrix(obj.size +1));  %Dimension in x−direction [m]
            Ly = str2double(inputMatrix(obj.size +2));  %Dimension in y−direction [m]
            h  = str2double(inputMatrix(obj.size +3));  %Reservoir thickness [m]
            Tres = str2double(inputMatrix(obj.temperature + 1));   %Res temperature [K]
            Reservoir = reservoir(Lx, Ly, h, Tres);
            phi = str2double(inputMatrix(obj.por + 1));
            if strcmp(inputMatrix(obj.perm - 1), 'INCLUDE')
                % File name
                file  = strcat('../Permeability/', char(inputMatrix(obj.perm +1)));
                % load the file in a vector
                field = load(file);
                % reshape it to specified size
                field = reshape(field(4:end),[field(1) field(2) field(3)]);
                % make it the size of the grid
                Kx = reshape(field(1:DiscretizationModel.ReservoirGrid.Nx,1:DiscretizationModel.ReservoirGrid.Ny, 1:DiscretizationModel.ReservoirGrid.Nz)*1e-15, DiscretizationModel.ReservoirGrid.N, 1);
                Ky = Kx;
                Kz = Kx;
            else
                value = str2double(inputMatrix(obj.perm +1));
                Kx = ones(DiscretizationModel.ReservoirGrid.N, 1)*value;
                Ky = ones(DiscretizationModel.ReservoirGrid.N, 1)*value;
                Kz = ones(DiscretizationModel.ReservoirGrid.N, 1)*value;
            end
            K = [Kx, Ky, Kz];
            Reservoir.AddPermeabilityPorosity(K, phi);
            ProductionSystem.AddReservoir(Reservoir);
            
            %% WELLS
            Wells = wells();
            Wells.NofInj = length(obj.inj);
            Wells.NofProd = length(obj.prod);
            n_phases = str2double(inputMatrix(obj.Comp_Type + 3));
            %Injectors
            for i=1:Wells.NofInj
                %i_init
                Well_Coord_Temp = strsplit(inputMatrix{obj.inj(i) + 1}, ' ');
                if sum( strcmp('NX' , Well_Coord_Temp) ) > 0
                    Well_Coord_Temp = strsplit(inputMatrix{obj.inj(i) + 1}, ' ');
                    if length(Well_Coord_Temp)>1
                        if Well_Coord_Temp{2}=='-',  i_init = DiscretizationModel.ReservoirGrid.Nx-1;
                        else,  error('For Injection Well Coordination, you can only use "NX" with minus sign "-"!');
                        end
                    else
                        i_init = DiscretizationModel.ReservoirGrid.Nx;   
                    end
                else
                    i_init = str2double(inputMatrix{obj.inj(i) + 1});
                end
                %i_final
                Well_Coord_Temp = strsplit(inputMatrix{obj.inj(i) + 2}, ' ');
                if sum( strcmp('NX' , Well_Coord_Temp) ) > 0
                    Well_Coord_Temp = strsplit(inputMatrix{obj.inj(i) + 2}, ' ');
                    if length(Well_Coord_Temp)>1
                        if Well_Coord_Temp{2}=='-',  i_final = DiscretizationModel.ReservoirGrid.Nx-1;
                        else,  error('For Injection Well Coordination, you can only use "NX" with minus sign "-"!');
                        end
                    else
                        i_final = DiscretizationModel.ReservoirGrid.Nx;   
                    end
                else
                    i_final = str2double(inputMatrix{obj.inj(i) + 2});
                end
                %j_init
                Well_Coord_Temp = strsplit(inputMatrix{obj.inj(i) + 3}, ' ');
                if sum( strcmp('NY' , Well_Coord_Temp) ) > 0
                    Well_Coord_Temp = strsplit(inputMatrix{obj.inj(i) + 3}, ' ');
                    if length(Well_Coord_Temp)>1
                        if Well_Coord_Temp{2}=='-',  j_init = DiscretizationModel.ReservoirGrid.Ny-1;
                        else,  error('For Injection Well Coordination, you can only use "NY" with minus sign "-"!');
                        end
                    else
                        j_init = DiscretizationModel.ReservoirGrid.Ny;   
                    end
                else
                    j_init = str2double(inputMatrix{obj.inj(i) + 3});
                end
                %j_final
                Well_Coord_Temp = strsplit(inputMatrix{obj.inj(i) + 4}, ' ');
                if sum( strcmp('NY' , Well_Coord_Temp) ) > 0
                    Well_Coord_Temp = strsplit(inputMatrix{obj.inj(i) + 4}, ' ');
                    if length(Well_Coord_Temp)>1
                        if Well_Coord_Temp{2}=='-',  j_final = DiscretizationModel.ReservoirGrid.Ny-1;
                        else,  error('For Injection Well Coordination, while you can only use "NY" with minus sign "-"!');
                        end
                    else
                        j_final = DiscretizationModel.ReservoirGrid.Ny;   
                    end
                else
                    j_final = str2double(inputMatrix{obj.inj(i) + 4});
                end
                %k_init
                Well_Coord_Temp = strsplit(inputMatrix{obj.inj(i) + 5}, ' ');
                if sum( strcmp('NZ' , Well_Coord_Temp) ) > 0
                    Well_Coord_Temp = strsplit(inputMatrix{obj.inj(i) + 5}, ' ');
                    if length(Well_Coord_Temp)>1
                        if Well_Coord_Temp{2}=='-',  k_init = DiscretizationModel.ReservoirGrid.Nz-1;
                        else,  error('For Injection Well Coordination, you can only use "NZ" with minus sign "-"!');
                        end
                    else
                        k_init = DiscretizationModel.ReservoirGrid.Nz;
                    end
                else
                    k_init = str2double(inputMatrix{obj.inj(i) + 5});
                end
                Well_Coord_Temp = strsplit(inputMatrix{obj.inj(i) + 6}, ' ');
                %k_final
                if sum( strcmp('NZ' , Well_Coord_Temp) ) > 0
                    Well_Coord_Temp = strsplit(inputMatrix{obj.inj(i) + 6}, ' ');
                    if length(Well_Coord_Temp)>1
                        if Well_Coord_Temp{2}=='-',  k_final = DiscretizationModel.ReservoirGrid.Nz-1;
                        else,  error('For Injection Well Coordination, while you can only use "NZ" with minus sign "-"!');
                        end
                    else
                        k_final = DiscretizationModel.ReservoirGrid.Nz;
                    end
                else
                    k_final = str2double(inputMatrix{obj.inj(i) + 6});
                end
                
                coord = [i_init, i_final; j_init, j_final; k_init, k_final];
                PI = 1000;
                pressure = str2double(inputMatrix(obj.inj(i) + 8));
                Injector = injector_pressure(PI, coord, pressure, Tres, n_phases);
                Wells.AddInjector(Injector);
            end
            
            %Producers
            for i=1:Wells.NofProd
                %i_init
                Well_Coord_Temp = strsplit(inputMatrix{obj.prod(i) + 1}, ' ');
                if sum( strcmp('NX' , Well_Coord_Temp) ) > 0
                    Well_Coord_Temp = strsplit(inputMatrix{obj.prod(i) + 1}, ' ');
                    if length(Well_Coord_Temp)>1
                        if Well_Coord_Temp{2}=='-',  i_init = DiscretizationModel.ReservoirGrid.Nx-1;
                        else,  error('For Production Well Coordination, you can only use "NX" with minus sign "-"!');
                        end
                    else
                        i_init = DiscretizationModel.ReservoirGrid.Nx;   
                    end
                else
                    i_init = str2double(inputMatrix{obj.prod(i) + 1});
                end
                %i_final
                Well_Coord_Temp = strsplit(inputMatrix{obj.prod(i) + 2}, ' ');
                if sum( strcmp('NX' , Well_Coord_Temp) ) > 0
                    Well_Coord_Temp = strsplit(inputMatrix{obj.prod(i) + 2}, ' ');
                    if length(Well_Coord_Temp)>1
                        if Well_Coord_Temp{2}=='-',  i_final = DiscretizationModel.ReservoirGrid.Nx-1;
                        else,  error('For Production Well Coordination, you can only use "NX" with minus sign "-"!');
                        end
                    else
                        i_final = DiscretizationModel.ReservoirGrid.Nx;   
                    end
                else
                    i_final = str2double(inputMatrix{obj.prod(i) + 2});
                end
                %j_init
                Well_Coord_Temp = strsplit(inputMatrix{obj.prod(i) + 3}, ' ');
                if sum( strcmp('NY' , Well_Coord_Temp) ) > 0
                    if length(Well_Coord_Temp)>1
                        if Well_Coord_Temp{2}=='-',  j_init = DiscretizationModel.ReservoirGrid.Ny-1;
                        else,  error('For Production Well Coordination, you can only use "NY" with minus sign "-"!');
                        end
                    else
                        j_init = DiscretizationModel.ReservoirGrid.Ny;   
                    end
                else
                    j_init = str2double(inputMatrix{obj.prod(i) + 3});
                end
                %j_final
                Well_Coord_Temp = strsplit(inputMatrix{obj.prod(i) + 4}, ' ');
                if sum( strcmp('NY' , Well_Coord_Temp) ) > 0
                    Well_Coord_Temp = strsplit(inputMatrix{obj.prod(i) + 4}, ' ');
                    if length(Well_Coord_Temp)>1
                        if Well_Coord_Temp{2}=='-',  j_final = DiscretizationModel.ReservoirGrid.Ny-1;
                        else,  error('For Production Well Coordination, while you can only use "NY" with minus sign "-"!');
                        end
                    else
                        j_final = DiscretizationModel.ReservoirGrid.Ny;   
                    end
                else
                    j_final = str2double(inputMatrix{obj.prod(i) + 4});
                end
                %k_init
                Well_Coord_Temp = strsplit(inputMatrix{obj.prod(i) + 5}, ' ');
                if sum( strcmp('NZ' , Well_Coord_Temp) ) > 0
                    Well_Coord_Temp = strsplit(inputMatrix{obj.prod(i) + 5}, ' ');
                    if length(Well_Coord_Temp)>1
                        if Well_Coord_Temp{2}=='-',  k_init = DiscretizationModel.ReservoirGrid.Nz-1;
                        else,  error('For Production Well Coordination, you can only use "NZ" with minus sign "-"!');
                        end
                    else
                        k_init = DiscretizationModel.ReservoirGrid.Nz;
                    end
                else
                    k_init = str2double(inputMatrix{obj.prod(i) + 5});
                end
                %k_final
                Well_Coord_Temp = strsplit(inputMatrix{obj.prod(i) + 6}, ' ');
                if sum( strcmp('NZ' , Well_Coord_Temp) ) > 0
                    Well_Coord_Temp = strsplit(inputMatrix{obj.prod(i) + 6}, ' ');
                    if length(Well_Coord_Temp)>1
                        if Well_Coord_Temp{2}=='-',  k_final = DiscretizationModel.ReservoirGrid.Nz-1;
                        else,  error('For prodection Well Coordination, while you can only use "NZ" with minus sign "-"!');
                        end
                    else
                        k_final = DiscretizationModel.ReservoirGrid.Nz;
                    end
                else
                    k_final = str2double(inputMatrix{obj.prod(i) + 6});
                end
                          
                coord = [i_init, i_final; j_init, j_final; k_init, k_final];
                PI = 1000;
                pressure = str2double(inputMatrix(obj.prod(i) + 8));
                Producer = producer_pressure(PI, coord, pressure);
                Wells.AddProducer(Producer);
            end
            ProductionSystem.AddWells(Wells);
            
            %% FRACTURE NETWORK
            if obj.Fractured == 1
                FracturesNetwork = fracture_system();
                FracturesNetwork.Active = 1;
                temp = strfind(FractureMatrix{1}, 'NUM_FRACS');
                index = find(~cellfun('isempty', temp));
                temp = strsplit(FractureMatrix{1}{index},' ');
                NrOfFrac = temp{end};
                FracturesNetwork.NumOfFrac = str2double( NrOfFrac );
                
                temp = strfind(FractureMatrix{1}, 'PROPERTIES');
                frac_index = find(~cellfun('isempty', temp));

                FracturesNetwork.Fractures = fracture();
                for f = 1 : FracturesNetwork.NumOfFrac
                    frac_info_split = strsplit(FractureMatrix{1}{frac_index(f)},' ');                   % Splitted data for each fracture
                    FracturesNetwork.Fractures(f).Length = str2double( frac_info_split{3} );            % Length of each fracture
                    FracturesNetwork.Fractures(f).Width = str2double( frac_info_split{4} );             % Width of each fracture
                    FracturesNetwork.Fractures(f).Thickness = str2double( frac_info_split{5} );         % Thickness of each fracture
                    
                    Porosity = str2double( frac_info_split{6} );                                        % Porosity  of each fracture 
                    Permeability = str2double( frac_info_split{7} );                                    % Permeability of each fracture
                    Kx = ones(DiscretizationModel.FracturesGrid.N(f), 1)*Permeability;
                    Ky = Kx;
                    Kz = Kx;
                    K = [Kx, Ky, Kz];
                    FracturesNetwork.Fractures(f).AddPermeabilityPorosity(K, Porosity);                 % Adding porosity and permeability to the fracture  
                end
                ProductionSystem.AddFractures(FracturesNetwork);
            end
        end
        function FluidModel = BuildFluidModel(obj, inputMatrix, ProductionSystem)
            n_phases = str2double(inputMatrix(obj.Comp_Type + 3));
            if n_phases == 1
                fluidmodel = 'SinglePhase';
                obj.CouplingType = 'SinglePhase';
            else
                fluidmodel = char(inputMatrix(obj.Comp_Type+1));
            end
            switch(fluidmodel)
                case('SinglePhase')
                    FluidModel = single_phase_fluid_model();
                    % Add phase
                    Phase = comp_phase();
                    %Gets all densities [kg/m^3]
                    Phase.rho0 = str2double(inputMatrix(obj.density + 2));
                    %Gets all viscosities [Pa sec]
                    Phase.mu = str2double(inputMatrix(obj.viscosity + 2));
                    % Compressibility
                    Phase.cf = str2double(inputMatrix(obj.compressibility + 2));
                    FluidModel.AddPhase(Phase, 1);
                    obj.Formulation = 'Immiscible';
                    if Phase.cf == 0
                        obj.incompressible = 1;
                    end
                case('Immiscible')
                    FluidModel = Immiscible_fluid_model(n_phases);
                    % Add phases
                    for i = 1:FluidModel.NofPhases
                        Phase = comp_phase();
                        %Gets all densities [kg/m^3]
                        Phase.rho0 = str2double(inputMatrix(obj.density + 2*i));
                        %Gets all viscosities [Pa sec]
                        Phase.mu = str2double(inputMatrix(obj.viscosity + 2*i));
                        % Compressibility
                        Phase.cf = str2double(inputMatrix(obj.compressibility + 2*i));
                        FluidModel.AddPhase(Phase, i);
                    end
                    obj.Formulation = 'Immiscible';
                case('BlackOil')
                    n_comp = n_phases;
                    % Black oil phases
                    FluidModel = BO_fluid_model(n_phases, n_comp);
                    % Gas
                    Gas = BO_gas_phase();
                    FluidModel.AddPhase(Gas, 1);
                    gas = BO_gas_component();
                    FluidModel.AddComponent(gas, 1);
                    % Oil
                    Oil = BO_oil_phase();
                    FluidModel.AddPhase(Oil, 2);
                    oil = BO_oil_component();
                    FluidModel.AddComponent(oil, 2);
                    if FluidModel.NofPhases == 3
                        Water = comp_phase();
                        Water.rho0 = 1000; % kg/m^3
                        Water.mu = 1e-3; % Pa s
                        Water.cf = 0;
                        FluidModel.AddPhase(Water, 3)
                        water = component();
                        FluidModel.AddComponent(water, 3);
                    end
                    FlashCalculator = Rachford_Rice_flash_calculator();
                    %FlashCalculator = Standard_flash_calculator();
                    FlashCalculator.KvaluesCalculator = BO_Kvalues_calculator();
                    FluidModel.FlashCalculator = FlashCalculator;
                case('Compositional')
                    n_comp = str2double(inputMatrix(obj.Comp_Type + 5));
                    FluidModel = Comp_fluid_model(n_phases, n_comp);
                    
                    %FlashCalculator = Standard_flash_calculator();
                    % Add phases
                    for i = 1:FluidModel.NofPhases
                        Phase = comp_phase();
                        %Gets all densities [kg/m^3]
                        Phase.rho0 = str2double(inputMatrix(obj.density + 2*i));
                        %Gets all viscosities [Pa sec]
                        Phase.mu = str2double(inputMatrix(obj.viscosity + 2*i));
                        % Compressibility
                        Phase.cf = str2double(inputMatrix(obj.compressibility + 2*i));
                        FluidModel.AddPhase(Phase, i);
                    end
                    % Add components
                    for i = 1:FluidModel.NofComp
                        Prop = str2double(strsplit(char(inputMatrix(obj.Comp_Prop + i * 2))));
                        comp = component();
                        comp.AddCompProperties(Prop);
                        FluidModel.AddComponent(comp, i);
                    end
                    FlashCalculator = Rachford_Rice_flash_calculator();
                    if length(Prop) > 2
                        FlashCalculator.KvaluesCalculator = Wilson_Kvalues_calculator(str2double(inputMatrix(obj.temperature + 1)));
                    else
                        FlashCalculator.KvaluesCalculator = Constant_Kvalues_calculator();
                    end
                    FluidModel.FlashCalculator = FlashCalculator;
            end
            
            %%  RelPerm model
            switch(char(inputMatrix(obj.relperm + 1)))
                case('Linear')
                    FluidModel.RelPermModel = relperm_model_linear();
                case('Quadratic')
                    FluidModel.RelPermModel = relperm_model_quadratic();
            end
            % Irriducible sat
            for i=1:FluidModel.NofPhases
                FluidModel.Phases(i).sr = str2double(inputMatrix(obj.relperm + 1 + 2*i));
            end
            %% Capillary pressure model
            switch (char(inputMatrix(obj.capillarity + 1)))
                case('JLeverett')
                    FluidModel.CapillaryModel = J_Function_model(ProductionSystem);
                    FluidModel.WettingPhaseIndex = str2double(inputMatrix(obj.capillarity + 2));
                case('Linear')
                    FluidModel.CapillaryModel = 'Not implemented';
                case('BrooksCorey')
                    FluidModel.CapillaryModel = 'Not implemented';
                case('Table')
                    FluidModel.CapillaryModel = 'Not implemented';
                otherwise
                    FluidModel.CapillaryModel = No_Pc_model();
            end
            
        end
        function Formulation = BuildFormulation(obj, Discretization, FluidModel,ProductionSystem)
            switch(obj.Formulation)
                case('Immiscible')
                    Formulation = Immiscible_formulation();
                    if strcmp(obj.ADM, 'active')
                        Discretization.OperatorsHandler.FullOperatorsAssembler = operators_assembler_Imm();
                    end
                case('Natural')
                    Formulation = NaturalVar_formulation(Discretization.ReservoirGrid.N, FluidModel.NofComp);
                    if strcmp(obj.ADM, 'active')
                        Discretization.OperatorsHandler.FullOperatorsAssembler = operators_assembler_comp();
                    end
                case('Molar')
                    Formulation = Overall_Composition_formulation(FluidModel.NofComp);
                    if strcmp(obj.ADM, 'active')
                        Discretization.OperatorsHandler.FullOperatorsAssembler = operators_assembler_Imm();
                    end
                case('OBL')
                    Formulation = OBL_formualtion();
                    Formulation.CreateTables();
            end
            Formulation.NofPhases = FluidModel.NofPhases;
            
            % Gravity model
            Formulation.GravityModel = gravity_model(Discretization, FluidModel.NofPhases, ProductionSystem.FracturesNetwork.NumOfFrac);
            switch (obj.Gravity)
                case('ON')
                    Formulation.GravityModel.g = 9.806;
                case('OFF')
                    Formulation.GravityModel.g = 0;
            end
        end
        function TimeDriver = BuildTimeDriver(obj, SettingsMatrix)
            TimeDriver = TimeLoop_Driver(obj.reports, obj.TotalTime, obj.MaxNumTimeSteps);
            %% Construct Coupling
            switch(obj.CouplingType)
                case('FIM')
                    %%%%FIM settings
                    % Build a different convergence cheker and a proper LS for ADM
                     NLSolver = NL_Solver();
                switch obj.ADM
                    case('inactive')
                        NLSolver.SystemBuilder = fim_system_builder();
                        switch (obj.Formulation)
                            case('Molar')
                                ConvergenceChecker = convergence_checker_FS_molar();
                            otherwise
                                ConvergenceChecker = convergence_checker_FS();
                        end
                        switch (obj.LinearSolver)
                            case ('gmres')
                                NLSolver.LinearSolver = linear_solver_iterative('gmres', 1e-6, 500);
                            case ('bicg')
                                NLSolver.LinearSolver = linear_solver_iterative('bicg', 1e-6, 500);
                            otherwise
                                NLSolver.LinearSolver = linear_solver();
                        end
                    case ('active')
                        NLSolver.SystemBuilder = fim_system_builder_ADM();
                        ConvergenceChecker = convergence_checker_ADM();
                        NLSolver.LinearSolver = linear_solver_ADM(obj.LinearSolver, 1e-6, 500); 
                end
                    NLSolver.MaxIter = str2double(SettingsMatrix(obj.coupling + 1));
                    ConvergenceChecker.Tol = str2double(SettingsMatrix(obj.coupling + 2));
                    switch (obj.Formulation)
                        case('Immiscible')
                            ConvergenceChecker.NormCalculator = norm_calculator_immiscible();
                        otherwise
                            ConvergenceChecker.NormCalculator = norm_calculator_comp();
                    end
                    NLSolver.AddConvergenceChecker(ConvergenceChecker);
                    % Build FIM Coupling strategy
                    Coupling = FIM_Strategy('FIM', NLSolver); 
                case('Sequential')
                    Coupling = Sequential_Strategy('Sequential');
                    Coupling.MaxIter = str2double(SettingsMatrix(obj.coupling + 1));
                    % pressuresolver = incompressible_pressure_solver();
                    pressuresolver = NL_Solver();
                    pressuresolver.MaxIter = 15;
                    ConvergenceChecker = convergence_checker_pressure();
                    ConvergenceChecker.Tol = 1e-6;
                    pressuresolver.AddConvergenceChecker(ConvergenceChecker);
                    pressuresolver.SystemBuilder = pressure_system_builder();
                    switch (obj.LinearSolver)
                            case ('gmres')
                                pressuresolver.LinearSolver = linear_solver_iterative('gmres', 1e-6, 500);
                            case ('bicg')
                                pressuresolver.LinearSolver = linear_solver_iterative('bicg', 1e-6, 500);
                            otherwise
                                pressuresolver.LinearSolver = linear_solver();
                    end
                    Coupling.AddPressureSolver(pressuresolver);
                    if (~isempty(obj.transport))
                        transportsolver = NL_Solver();
                        transportsolver.MaxIter = str2double(SettingsMatrix(obj.transport + 3));
                        ConvergenceChecker = convergence_checker_transport();
                        ConvergenceChecker.Tol = str2double(SettingsMatrix(obj.transport + 2));
                        transportsolver.AddConvergenceChecker(ConvergenceChecker);
                        transportsolver.SystemBuilder = transport_system_builder();
                        Coupling.ConvergenceChecker = convergence_checker_outer();
                        Coupling.ConvergenceChecker.Tol = str2double(SettingsMatrix(obj.coupling + 2));
                        switch (obj.LinearSolver)
                            case ('gmres')
                                transportsolver.LinearSolver = linear_solver_iterative('gmres', 1e-6, 500);
                            case ('bicg')
                                transportsolver.LinearSolver = linear_solver_iterative('bicg', 1e-6, 500);
                            otherwise
                                transportsolver.LinearSolver = linear_solver();
                        end
                    else
                        transportsolver = explicit_transport_solver();
                        transportsolver.SystemBuilder = explicit_transport_system_builder();
                        Coupling.ConvergenceChecker = convergence_checker_impes();
                    end
                    Coupling.AddTransportSolver(transportsolver);
                case('SinglePhase')
                    Coupling = SinglePhase_Strategy('SinglePhase');
                    pressuresolver = NL_Solver();
                    pressuresolver.MaxIter = 15;
                    ConvergenceChecker = convergence_checker_pressure();
                    ConvergenceChecker.Tol = 1e-6;
                    pressuresolver.AddConvergenceChecker(ConvergenceChecker);
                    pressuresolver.SystemBuilder = pressure_system_builder();    
                    pressuresolver.LinearSolver = linear_solver();
                    Coupling.AddPressureSolver(pressuresolver);
                    if obj.incompressible
                        Coupling.Incompressible = 1;
                    end
            end
            Coupling.TimeStepSelector = timestep_selector(str2double(SettingsMatrix(obj.coupling + 3)));
            TimeDriver.AddCouplingStrategy(Coupling);
            switch(obj.StopCriterion)
                case('MAX TIME')
                    end_of_sim_eval = end_of_sim_evaluator(obj.TotalTime, obj.MaxNumTimeSteps);
                case('COMPONENT CUT')
                    end_of_sim_eval = end_of_sim_evaluator_gascut(obj.TotalTime, obj.MaxNumTimeSteps, 0.2);
                case('PV INJECTED')
                    end_of_sim_eval = end_of_sim_evaluator_PVInjected(obj.TotalTime, obj.MaxNumTimeSteps, 3);
            end
            TimeDriver.AddEndOfSimEvaluator(end_of_sim_eval);
        end
        function Summary = BuildSummary(obj, simulation)
            %%%%%%%%%%%%%%% BuildObjects for OUTPUT%%%%%%%%%
            switch(obj.CouplingType)
                case('FIM')
                    CouplingStats = FIM_Stats(obj.MaxNumTimeSteps);
                case('Sequential')
                    CouplingStats = Sequential_Stats(obj.MaxNumTimeSteps);
                case('SinglePhase')
                    CouplingStats = SinglePhase_Stats(obj.MaxNumTimeSteps);
            end
            wellsData = wells_data(obj.MaxNumTimeSteps, simulation.FluidModel.NofPhases, simulation.FluidModel.NofComp, simulation.ProductionSystem.Wells);
            switch (obj.ADM)
                case('inactive')
                    Summary = Run_Summary(obj.MaxNumTimeSteps, CouplingStats, wellsData);
                case('active')
                    Summary = Run_Summary_ADM(obj.MaxNumTimeSteps, CouplingStats, wellsData, simulation.DiscretizationModel.maxLevel(1)); % Only reservoir for now
            end
        end
        function Writer = BuildWriter(obj, InputDirectory, simulation)
            % Build Plotter
            switch(obj.plotting)
                case('Matlab')
                    if simulation.DiscretizationModel.ReservoirGrid.Nx == 1 || simulation.DiscretizationModel.ReservoirGrid.Ny == 1
                        plotter = Matlab_Plotter_1D();
                    else
                        plotter = Matlab_Plotter_2D();
                    end
                case('VTK')
                    plotter = VTK_Plotter(InputDirectory, obj.ProblemName);
                otherwise
                    warning('WARNING: NO valid Plotter was selected. Results will not be plotted.');
                    plotter = no_Plotter();
            end
            
            switch(obj.ADM)
                case ('inactive')
                    Writer = output_writer_FS(InputDirectory, obj.ProblemName, simulation.ProductionSystem.Wells.NofInj, simulation.ProductionSystem.Wells.NofProd, simulation.Summary.CouplingStats.NTimers, simulation.Summary.CouplingStats.NStats, simulation.FluidModel.NofComp);
                case ('active')
                    Writer = output_writer_adm(InputDirectory, obj.ProblemName, simulation.ProductionSystem.Wells.NofInj, simulation.ProductionSystem.Wells.NofProd, simulation.Summary.CouplingStats.NTimers, simulation.Summary.CouplingStats.NStats, simulation.FluidModel.NofComp);
            end
            Writer.AddPlotter(plotter);
        end
        function DefineProperties(obj, ProductionSystem, FluidModel, DiscretizationModel)
            switch(obj.Fractured)
                case(0)
                    ProductionSystem.Reservoir.State.AddProperties(FluidModel, DiscretizationModel.ReservoirGrid.N);
                    ProductionSystem.Reservoir.State_old.AddProperties(FluidModel, DiscretizationModel.ReservoirGrid.N);
                case(1)
                    ProductionSystem.Reservoir.State.AddProperties(FluidModel, DiscretizationModel.ReservoirGrid.N);
                    ProductionSystem.Reservoir.State_old.AddProperties(FluidModel, DiscretizationModel.ReservoirGrid.N);
                    for f=1:ProductionSystem.FracturesNetwork.NumOfFrac
                        ProductionSystem.FracturesNetwork.Fractures(f).State.AddProperties(FluidModel, DiscretizationModel.FracturesGrid.N(f));
                        ProductionSystem.FracturesNetwork.Fractures(f).State_old.AddProperties(FluidModel, DiscretizationModel.FracturesGrid.N(f));
                    end
            end
        end
    end
end