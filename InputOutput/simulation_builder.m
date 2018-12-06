 % Builder Builds all objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef simulation_builder < handle
    properties
        SimulationInput
        SimulatorSettings
        incompressible = 0;
        NofEq
    end
    methods
        function simulation = BuildSimulation(obj, FractureMatrix)
            simulation = Reservoir_Simulation();
            simulation.DiscretizationModel = obj.BuildDiscretization(FractureMatrix);
            simulation.ProductionSystem = obj.BuildProductionSystem(FractureMatrix,simulation.DiscretizationModel.FracturesGrid);            
            simulation.FluidModel = obj.BuildFluidModel();
            simulation.Formulation = obj.BuildFormulation();
            simulation.TimeDriver = obj.BuildTimeDriver();
            simulation.Summary = obj.BuildSummary(simulation);
            
            % Add gravity model
            simulation.Formulation.GravityModel = gravity_model(simulation.DiscretizationModel, simulation.FluidModel.NofPhases, simulation.ProductionSystem.FracturesNetwork.NumOfFrac);
            switch (obj.SimulationInput.FluidProperties.Gravity)
                case('ON')
                    simulation.Formulation.GravityModel.g = 9.806;
                case('OFF')
                    simulation.Formulation.GravityModel.g = 0;
            end
            
            % Define Properties
            obj.DefineProperties(simulation.ProductionSystem, simulation.FluidModel, simulation.DiscretizationModel);
            
            %% Define Initialization procedure
            N = simulation.DiscretizationModel.ReservoirGrid.N;
            VarValues = ones(N, length(obj.SimulationInput.Init));
            for i=1:length(obj.SimulationInput.Init)
                VarValues(:, i) = VarValues(:, i) * obj.SimulationInput.Init(i);
            end
            switch(simulation.FluidModel.name)
                case('SinglePhase') 
                    VarNames = {'P_1', 'S_1'};
                    VarValues(:, 2) = 1;
                    simulation.Initializer = initializer_singlephase(VarNames, VarValues);
                case('Immiscible')
                    VarNames = {'P_2', 'S_1', 'S_2'};
%                     % Perturb initial S
%                     nx = simulation.DiscretizationModel.ReservoirGrid.Nx;
%                     rng(0);
%                     perturbations = rand(20, 1);
%                     np = nx / 20;
%                     for i=1:20
%                         newval((i-1)*np + 1:i*np) = perturbations(i);
%                     end
%                     index = 1:nx:N;
%                     VarValues(index, 2) = newval;
%                     VarValues(:, 3) = 1 - VarValues(:, 2);
                    simulation.Initializer = initializer(VarNames, VarValues);
                otherwise
                    VarNames = {'P_2', 'z_1', 'z_2'};
                    simulation.Initializer = initializer_hydrostatic(VarNames, VarValues);
            end
        end
        function Discretization = BuildDiscretization(obj, FractureMatrix)
            %% 1. Create fine-scale grids
            % 1a. Reservoir Grid
            ReservoirGrid = cartesian_grid(obj.SimulationInput.ReservoirProperties.Grid.N);
            Nm = ReservoirGrid.N;
            % 1b. Fractures Grid
            if obj.SimulationInput.FracturesProperties.Fractured
                NrOfFrac = obj.SimulationInput.FracturesProperties.NrOfFrac;
                fprintf('Extracting data from %02d fractures ...\n', NrOfFrac);
                FracturesGrid = fractures_grid(NrOfFrac);
                temp = strfind(FractureMatrix, 'PROPERTIES');
                frac_index = find(~cellfun('isempty', temp));
                
                Nf = zeros(NrOfFrac,1);
                for f = 1 : NrOfFrac
                    % Creat cartesian grid in each fracture
                    frac_info_split = strsplit(FractureMatrix{frac_index(f)},' ');
                    grid_temp = strsplit(frac_info_split{8}, 'x');
                    nx = str2double(grid_temp{1});
                    ny = str2double(grid_temp{2});
                    nz = 1;
                    FractureGrid = cartesian_grid([nx;ny;nz]);
                    
                    % Add fracture grid coordinates (it's for plotting purposes) 
                    temp = strfind(FractureMatrix, 'GRID_COORDS_X');
                    frac_grid_coords_x = find(~cellfun('isempty', temp));
                    temp = strfind(FractureMatrix, 'GRID_COORDS_Y');
                    frac_grid_coords_y = find(~cellfun('isempty', temp));
                    temp = strfind(FractureMatrix, 'GRID_COORDS_Z');
                    frac_grid_coords_z = find(~cellfun('isempty', temp));
                    frac_grid_coords_x_split = strsplit(FractureMatrix{frac_grid_coords_x(f)},' ');
                    frac_grid_coords_y_split = strsplit(FractureMatrix{frac_grid_coords_y(f)},' ');
                    frac_grid_coords_z_split = strsplit(FractureMatrix{frac_grid_coords_z(f)},' ');
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
                temp = strfind(FractureMatrix, 'FRACCELL');
                frac_cell_index = find(~cellfun('isempty', temp));
                
                temp = strfind(FractureMatrix, 'ROCK_CONN');
                frac_rockConn_index = find(~cellfun('isempty', temp));
                
                temp = strfind(FractureMatrix, 'FRAC_CONN');
                frac_fracConn_index = find(~cellfun('isempty', temp));
                
                %n_phases = str2double(inputMatrix(obj.Comp_Type + 3)); % Number of phases (useful to define size of some objects)
                n_phases = obj.SimulationInput.FluidProperties.NofPhases;
                fprintf('---> Fracture ');
                for f = 1 : NrOfFrac
                    if (f>1),  fprintf(repmat('\b', 1, 5+27));  end
                    fprintf('%02d/%02d',f,NrOfFrac);
                    % looping over all global fracture cells
					fprintf(' ---> Grid cell ');
                    for If = 1:Nf(f)
						if (If>1),  fprintf(repmat('\b', 1, 11));  end
						fprintf('%05d/%05d',If,Nf(f));
                        fracCell_info_split = strsplit(FractureMatrix{frac_cell_index(sum(Nf(1:f-1))+If)},{' ','	'});
                        
                        % frac-marix conn
                        temp = frac_rockConn_index - frac_cell_index(sum(Nf(1:f-1))+If); temp(temp<0) = max(temp) +1;
                        [~ , frac_rockConn_index_start] = min(temp);
                        for Im = 1:str2double(fracCell_info_split{3})
                            frac_rockConn_info_split = strsplit(FractureMatrix{frac_rockConn_index(frac_rockConn_index_start+Im-1)},{' ','	'});
                            CrossConnections(sum(Nf(1:f-1))+If,1).Cells(Im,1) = str2double(frac_rockConn_info_split{2})+1;
                            CrossConnections(sum(Nf(1:f-1))+If,1).ConnIndex(Im,1) = str2double(frac_rockConn_info_split{3});
                        end
                        
                        % frac-frac conn
                        temp = frac_fracConn_index - frac_cell_index(sum(Nf(1:f-1))+If); temp(temp<0) = max(temp)+1;
                        [~ , frac_fracConn_index_start] = min(temp);
                        
                        Counter = 1;
                        for Ig = 1:str2double(fracCell_info_split{5})
                            frac_fracConn_info_split = strsplit(FractureMatrix{frac_fracConn_index(frac_fracConn_index_start+Ig-1)},{' ','	'});
                            If_Other_Global = Nm + sum( Nf( 1 : str2double(frac_fracConn_info_split{2}) +1-1 ) ) + str2double(frac_fracConn_info_split{3})+1;
                            if If_Other_Global > sum(Nf(1:f-1))+If + Nm
                                CrossConnections(sum(Nf(1:f-1))+If,1).Cells(Im+Counter,1) = If_Other_Global;
                                CrossConnections(sum(Nf(1:f-1))+If,1).ConnIndex(Im+Counter,1) = str2double(frac_fracConn_info_split{4});
                                Counter = Counter + 1;
                            end
                        end
                        CrossConnections(sum(Nf(1:f-1))+If,1).UpWind = zeros(length(CrossConnections(sum(Nf(1:f-1))+If,1).Cells), n_phases);
                        CrossConnections(sum(Nf(1:f-1))+If,1).U_Geo = zeros(length(CrossConnections(sum(Nf(1:f-1))+If,1).Cells), n_phases);
                        CrossConnections(sum(Nf(1:f-1))+If,1).T_Geo = zeros(size(CrossConnections(sum(Nf(1:f-1))+If,1).ConnIndex));
                    end
                end
                fprintf(' ---> Completed');
				fprintf('\n');
            end
            
            %% 2. Define your discretization Model (choose between FS and ADM)
            switch (obj.SimulatorSettings.DiscretizationModel)
                case('ADM')
                    % ADM discretization model (dynamic gridding)
                    ADMSettings = obj.SimulatorSettings.ADMSettings;
                    % Create the operatorshandler
                    operatorshandler = operators_handler_adm(ADMSettings.Coarsening(1,:,:));
                    % a.1 Pressure prolongation builder
                    switch (ADMSettings.PInterpolator)
                        case ('Constant')
                            prolongationbuilder = prolongation_builder_constant(ADMSettings.maxLevel(1));
                        otherwise
                            prolongationbuilder = prolongation_builder_MSPressure(ADMSettings.maxLevel(1), ADMSettings.Coarsening(:,:,1));
                            if ~obj.SimulationInput.FracturesProperties.Fractured
                                prolongationbuilder.BFUpdater = bf_updater_ms();
                            else
                                prolongationbuilder.BFUpdater = bf_updater_FAMS();
                                prolongationbuilder.BFUpdater.BFtype = ADMSettings.BFtype;
                            end
                            if strcmp(ADMSettings.PInterpolator, 'Homogeneous')
                                prolongationbuilder.BFUpdater.MaxContrast = 1;
                            else
                                prolongationbuilder.BFUpdater.MaxContrast = 10^-2;
                            end
                    end
                    
                    operatorshandler.AddProlongationBuilder(prolongationbuilder, 1);
                    % a.2 Hyperbolic variables operators builder
                    n_phases = obj.SimulationInput.FluidProperties.NofPhases;
                    for i = 2:n_phases
                        switch(ADMSettings.HInterpolator)
                            case('Constant')
                                prolongationbuilder = prolongation_builder_constant(ADMSettings.maxLevel(1));
                            case('MS')
                                prolongationbuilder = prolongation_builder_MSHyperbolic(ADMSettings.maxLevel(1));
                        end
                        operatorshandler.AddProlongationBuilder(prolongationbuilder, i);
                    end
                    
                    % b. Grid selection criterion (time\space based)
                    switch (ADMSettings.GridSelCriterion)
                        case('dfdx')
                            gridselector = adm_grid_selector_delta(ADMSettings.tol, ADMSettings.key);
                        case('dfdt')
                            gridselector = adm_grid_selector_time(ADMSettings.tol, ADMSettings.key, ReservoirGrid.N, ADMSettings.maxLevel(1));
                        case('residual')
                            gridselector = adm_grid_selector_residual(ADMSettings.tol);
                    end
                    Discretization = ADM_Discretization_model(ADMSettings.maxLevel, ADMSettings.Coarsening);
                    Discretization.AddADMGridSelector(gridselector);
                    Discretization.AddOperatorsHandler(operatorshandler);
                case ('MMs')
                    % Multilevel multiscale discretization model (static multiscale)
                    MMsSettings = obj.SimulatorSettings.MMsSettings;
                    % Create the operatorshandler
                    operatorshandler = operators_handler_MMs(MMsSettings.Coarsening(1,:,:));
                    prolongationbuilder = prolongation_builder_MSPressure(MMsSettings.maxLevel(1), MMsSettings.Coarsening(:,:,1) );
                    if ~obj.SimulationInput.FracturesProperties.Fractured
                        prolongationbuilder.BFUpdater = bf_updater_ms();
                    else
                        prolongationbuilder.BFUpdater = bf_updater_FAMS();
                        if strcmp(MMsSettings.BFtype , 'COUPLED')
                            prolongationbuilder.BFUpdater.BFtype = 'COUPLED';
                        else
                            prolongationbuilder.BFUpdater.BFtype = 'DECOUPLED';
                        end
                    end
                    % Reduce contrast for BF computation to remove peaks
                    prolongationbuilder.BFUpdater.MaxContrast = 10^-2;
                    if MMsSettings.CorrectionFunctions
                        prolongationbuilder.BFUpdater.CorrectionFunctions = true;
                    end
                    
                    % Add prolongation builder
                    operatorshandler.AddProlongationBuilder(prolongationbuilder, 1);
                    % Static Multiscale for flow solver
                    Discretization = Multiscale_Discretization_model(MMsSettings.maxLevel, MMsSettings.Coarsening);
                    Discretization.AddOperatorsHandler(operatorshandler);
                case('FS')
                    % Fine-scale discretization model
                    Discretization = FS_Discretization_model();
            end
            
            %% 3. Add Grids to the Discretization Model
            Discretization.AddReservoirGrid(ReservoirGrid);
            if obj.SimulationInput.FracturesProperties.Fractured
                Discretization.AddFracturesGrid(FracturesGrid);
                Discretization.AddCrossConnections(CrossConnections);
            end
        end
        function ProductionSystem = BuildProductionSystem (obj, FractureMatrix, FracturesGrid)
            ProductionSystem = Production_System();
            %% RESERVOIR
            Lx = obj.SimulationInput.ReservoirProperties.size(1);       %Dimension in x-direction [m]
            Ly = obj.SimulationInput.ReservoirProperties.size(2);       %Dimension in y-direction [m]
            h  = obj.SimulationInput.ReservoirProperties.size(3);       %Reservoir thickness (z-direction) [m]
            Tres = obj.SimulationInput.ReservoirProperties.Temperature; %Res temperature [K]
            Reservoir = reservoir(Lx, Ly, h, Tres);
            phi = obj.SimulationInput.ReservoirProperties.phi;
            nx = obj.SimulationInput.ReservoirProperties.Grid.N(1); 
            ny = obj.SimulationInput.ReservoirProperties.Grid.N(2);
            nz = obj.SimulationInput.ReservoirProperties.Grid.N(3);
            K = ones(nx*ny*nz, 3);
            for i=1:3
                if obj.SimulationInput.ReservoirProperties.PermInclude(i)
                    % load the file in a vector
                    if i==1
                        fprintf('\n---> Reading permeability file #%d ...',i);
                        field = load(obj.SimulationInput.ReservoirProperties.PermFile{i});
                        fprintf(' ---> Completed.\n');
                    else
                        % loading the permeability file only if different
                        % than previous one
                        if ~strcmp(obj.SimulationInput.ReservoirProperties.PermFile{i},obj.SimulationInput.ReservoirProperties.PermFile{i-1})
                            fprintf('---> Reading permeability file #%d ...',i);
                            field = load(obj.SimulationInput.ReservoirProperties.PermFile{i});
                            fprintf(' ---> Completed.\n');
                        end
                    end
                    % reshape it to specified size
                    field1 = reshape(field(4:end,1),[field(1,1) field(2,1) field(3,1)]);
                    % make it the size of the grid
                    %K(:,i) = reshape(field1(1:nx, 1:ny, 1:nz)*1e-15, nx*ny*nz, 1);
                    % In case the data is in logarithmic scale
                    K(:,i) = reshape(field1(1:nx, 1:ny, 1:nz)*1e-15, nx*ny*nz, 1);
                else
                    value = obj.SimulationInput.ReservoirProperties.Perm(i);
                    K(:, i)= K(:,i) * value;
                end
            end
            Reservoir.AddPermeabilityPorosity(K, phi);
            switch obj.SimulatorSettings.DiscretizationModel
                case('ADM')
                    % This is for DLGR type ADM: it reads coarse permeabilities
                    if obj.SimulatorSettings.ADMSettings.DLGR
                        K_coarse = cell(obj.SimulatorSettings.ADMSettings.maxLevel + 1, 1);
                        K_coarse{1} = K;
                        for l=2:obj.SimulatorSettings.ADMSettings.maxLevel + 1
                            for d=1:2
                                % load the file in a vector
                                field = load(obj.SimulationInput.ReservoirProperties.CoarsePermFile{l-1,d});
                                % reshape it to specified size
                                k = field(4:end)*1e-15;
                                K_coarse{l}(:, d) = k;
                            end
                            K_coarse{l}(:, 3) = k;
                        end
                        % Save them in ProductionSystem.
                        Reservoir.AddCoarsePermeability(K_coarse); % this function you have to create it
                    end
                otherwise
            end
            % Add reservoir to production system            
            ProductionSystem.AddReservoir(Reservoir);
            
            %% WELLS
            Wells = wells();
            Wells.NofInj = obj.SimulationInput.WellsInfo.NofInj;
            Wells.NofProd = obj.SimulationInput.WellsInfo.NofProd;
            n_phases = obj.SimulationInput.FluidProperties.NofPhases;
            
            % Useful to compute PI
            dx = Lx/nx;
            dy = Ly/ny;
            dz = h /nz;
            %Injectors
            for i=1:Wells.NofInj
                switch(obj.SimulationInput.WellsInfo.Inj(i).PI.type)
                    case ('Dirichlet')
                        PI = dy*dz/(dx/2);
                    case ('PI')
                        PI = obj.SimulationInput.WellsInfo.Inj(i).PI.value;
                    case ('Radius')
                        radius = obj.SimulationInput.WellsInfo.Inj(i).PI.value;
                        error('DARSim2 error: Radius calculation of PI is not implemented for now')
                end
                coord = obj.SimulationInput.WellsInfo.Inj(i).Coord;
                switch (obj.SimulationInput.WellsInfo.Inj(i).Constraint.name)
                    case('pressure')
                        pressure = obj.SimulationInput.WellsInfo.Inj(i).Constraint.value;
                        temperature = obj.SimulationInput.WellsInfo.Inj(i).Temperature;
                        Injector = injector_pressure(PI, coord, pressure, temperature, n_phases);
                    case('rate')
                        rate = obj.SimulationInput.WellsInfo.Inj(i).Constraint.value;
                        p_init = obj.SimulationInput.Init(1);
                        rate = rate * Reservoir.TotalPV / (3600 * 24); % convert pv/day to m^3/s
                        Injector = injector_rate(PI, coord, rate, p_init, Tres, n_phases);
                end
                Wells.AddInjector(Injector);
            end
            
            %Producers
            for i=1:Wells.NofProd
                switch(obj.SimulationInput.WellsInfo.Prod(i).PI.type)
                    case ('Dirichlet')
                        PI = dy*dz/(dx/2);
                    case ('PI')
                        PI = obj.SimulationInput.WellsInfo.Prod(i).PI.value;
                    case ('Radius')
                        radius = obj.SimulationInput.WellsInfo.Prod(i).PI.value;
                        error('DARSim2 error: Radius calculation of PI is not implemented for now')
                end
                coord = obj.SimulationInput.WellsInfo.Prod(i).Coord;
                switch (obj.SimulationInput.WellsInfo.Prod(i).Constraint.name)
                    case('pressure')
                        pressure = obj.SimulationInput.WellsInfo.Prod(i).Constraint.value;
                        Producer = producer_pressure(PI, coord, pressure);
                    case('rate')
                        rate = obj.SimulationInput.WellsInfo.Prod(i).Constraint.value;
                        rate = rate * Reservoir.TotalPV / (3600 * 24); % convert pv/day to m^3/s
                        Producer = producer_rate(PI, coord, rate);
                end
                Wells.AddProducer(Producer);
            end
            ProductionSystem.AddWells(Wells);
            
            %% FRACTURE NETWORK
            if obj.SimulationInput.FracturesProperties.Fractured
                FracturesNetwork = fracture_system();
                FracturesNetwork.Active = 1;
                temp = strfind(FractureMatrix, 'NUM_FRACS');
                index = find(~cellfun('isempty', temp));
                temp = strsplit(FractureMatrix{index},' ');
                NrOfFrac = temp{end};
                FracturesNetwork.NumOfFrac = str2double( NrOfFrac );
                
                temp = strfind(FractureMatrix, 'PROPERTIES');
                frac_index = find(~cellfun('isempty', temp));

                FracturesNetwork.Fractures = fracture();
                for f = 1 : FracturesNetwork.NumOfFrac
                    frac_info_split = strsplit(FractureMatrix{frac_index(f)},' ');                      % Splitted data for each fracture
                    FracturesNetwork.Fractures(f).Length = str2double( frac_info_split{3} );            % Length of each fracture
                    FracturesNetwork.Fractures(f).Width = str2double( frac_info_split{4} );             % Width of each fracture
                    FracturesNetwork.Fractures(f).Thickness = str2double( frac_info_split{5} );         % Thickness of each fracture
                    
                    Porosity = str2double( frac_info_split{6} );                                        % Porosity  of each fracture 
                    Permeability = str2double( frac_info_split{7} );                                    % Permeability of each fracture
                    Kx = ones(FracturesGrid.N(f), 1)*Permeability;
                    Ky = Kx;
                    Kz = Kx;
                    K = [Kx, Ky, Kz];
                    FracturesNetwork.Fractures(f).AddPermeabilityPorosity(K, Porosity);                 % Adding porosity and permeability to the fracture  
                end
                ProductionSystem.AddFractures(FracturesNetwork);
            end
        end
        function FluidModel = BuildFluidModel(obj)
            n_phases = obj.SimulationInput.FluidProperties.NofPhases;
            switch(obj.SimulationInput.FluidProperties.FluidModel)
                case('SinglePhase')
                    FluidModel = single_phase_fluid_model();
                    % Add phase
                    Phase = comp_phase();
                    %Gets all densities [kg/m^3]
                    Phase.rho0 = obj.SimulationInput.FluidProperties.Density;
                    %Gets all viscosities [Pa sec]
                    Phase.mu = obj.SimulationInput.FluidProperties.mu;
                    % Compressibility
                    Phase.cf = obj.SimulationInput.FluidProperties.Comp;
                    FluidModel.AddPhase(Phase, 1);
                    if Phase.cf == 0
                        obj.incompressible = 1;
                    end
                    obj.SimulatorSettings.CouplingType = 'SinglePhase';
                case('Immiscible')
                    FluidModel = Immiscible_fluid_model(n_phases);
                    % Add phases
                    for i = 1:FluidModel.NofPhases
                        Phase = comp_phase();
                        %Gets all densities [kg/m^3]
                        Phase.rho0 = obj.SimulationInput.FluidProperties.Density(i);
                        %Gets all viscosities [Pa sec]
                        Phase.mu = obj.SimulationInput.FluidProperties.mu(i);
                        % Compressibility
                        Phase.cf = obj.SimulationInput.FluidProperties.Comp(i);
                        FluidModel.AddPhase(Phase, i);
                    end
                    obj.SimulatorSettings.Formulation = 'Immiscible';
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
                    n_comp = obj.SimulationInput.FluidProperties.NofComp;
                    FluidModel = Comp_fluid_model(n_phases, n_comp);
                    %FlashCalculator = Standard_flash_calculator();
                    % Add phases
                    for i = 1:FluidModel.NofPhases
                        Phase = comp_phase();
                        %Gets all densities [kg/m^3]
                        Phase.rho0 = obj.SimulationInput.FluidProperties.Density(i);
                        %Gets all viscosities [Pa sec]
                        Phase.mu = obj.SimulationInput.FluidProperties.mu(i);
                        % Compressibility
                        Phase.cf = obj.SimulationInput.FluidProperties.Comp(i);
                        FluidModel.AddPhase(Phase, i);
                    end
                    % Add components
                    for i = 1:FluidModel.NofComponents
                        Prop = obj.SimulationInput.FluidProperties.ComponentProperties(i,:);
                        comp = component();
                        comp.AddCompProperties(Prop);
                        FluidModel.AddComponent(comp, i);
                    end
                    FlashCalculator = Rachford_Rice_flash_calculator();
                    if length(Prop) > 2
                        FlashCalculator.KvaluesCalculator = Wilson_Kvalues_calculator(obj.SimulationInput.ReservoirProperties.Temperature);
                    else
                        FlashCalculator.KvaluesCalculator = Constant_Kvalues_calculator();
                    end
                    FluidModel.FlashCalculator = FlashCalculator;
                case('Geothermal')
                    % build the geothermal fluid model
                    FluidModel = Geothermal_fluid_model();
                    Phase = therm_comp_phase();
                    FluidModel.AddPhase(Phase, 1);
                    obj.SimulatorSettings.CouplingType = 'FIM';
            end
            
            %%  RelPerm model
            switch(obj.SimulationInput.FluidProperties.RelPerm.name)
                case('Linear')
                    FluidModel.RelPermModel = relperm_model_linear();
                case('Quadratic')
                    FluidModel.RelPermModel = relperm_model_quadratic();
                case('Foam')
                    FluidModel.RelPermModel = relperm_model_foam();
                case('Corey')
                    FluidModel.RelPermModel = relperm_model_brookscorey();
            end
            % Irriducible sat
            for i=1:FluidModel.NofPhases
                FluidModel.Phases(i).sr = obj.SimulationInput.FluidProperties.RelPerm.s_irr(i);
            end
            %% Capillary pressure model
            switch (obj.SimulationInput.FluidProperties.Capillarity.name)
                case('JLeverett')
                    FluidModel.CapillaryModel = J_Function_model();
                    FluidModel.WettingPhaseIndex = obj.SimulationInput.FluidProperties.Capillarity.wetting;
                case('Linear')
                    FluidModel.CapillaryModel = 'Not implemented';
                case('Corey')
                    FluidModel.CapillaryModel = 'Not implemented';
                case('Table')
                    FluidModel.CapillaryModel = 'Not implemented';
                otherwise
                    FluidModel.CapillaryModel = No_Pc_model();
            end
            
        end
        function Formulation = BuildFormulation(obj)
            switch(obj.SimulatorSettings.Formulation)
                case('Immiscible')
                    Formulation = Immiscible_formulation();
                    obj.NofEq = obj.SimulationInput.FluidProperties.NofPhases;
                case('Natural')
                    Formulation = NaturalVar_formulation(sum(obj.SimulationInput.ReservoirProperties.Grid.N), obj.SimulationInput.FluidProperties.NofComponents);
                    obj.NofEq = obj.SimulationInput.FluidProperties.NofPhases + obj.SimulationInput.FluidProperties.NofComponents;
                case('Molar')
                    Formulation = Overall_Composition_formulation(obj.SimulationInput.FluidProperties.NofComponents);
                    obj.NofEq = obj.SimulationInput.FluidProperties.NofPhases;
                case('OBL')
                    Formulation = OBL_formualtion();
                    Formulation.CreateTables();
                case('Thermal')
                    Formulation = Thermal_formulation();
                    obj.NofEq = obj.SimulationInput.FluidProperties.NofPhases + 1;
            end
            Formulation.NofPhases = obj.SimulationInput.FluidProperties.NofPhases;
        end
        function TimeDriver = BuildTimeDriver(obj)
            TimeDriver = TimeLoop_Driver(obj.SimulatorSettings.reports, obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps);
            % Construct Coupling
            switch(obj.SimulatorSettings.CouplingType)
                case('FIM')
                    % FIM coupling
                    %%%%FIM settings
                    NLSolver = NL_Solver();
                    switch obj.SimulatorSettings.DiscretizationModel
                        case ('ADM')
                            % Build a different convergence cheker and a proper LS for ADM
                            NLSolver.SystemBuilder = fim_system_builder_ADM();
                            ConvergenceChecker = convergence_checker_ADM();
                            ConvergenceChecker.OperatorsAssembler = operators_assembler_fim(obj.NofEq);
                            NLSolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            NLSolver.LinearSolver.OperatorsAssembler = operators_assembler_fim(obj.NofEq);
                            if obj.SimulatorSettings.ADMSettings.DLGR
                                % it will change perm during ADM
                                % simulaiton to use upscaled ones
                                NLSolver.LinearSolver.DLGR = 1;
                            end
                        otherwise
                            NLSolver.SystemBuilder = fim_system_builder();
                            switch (obj.SimulatorSettings.Formulation)
                                case('Molar')
                                    ConvergenceChecker = convergence_checker_FS_molar();
                                otherwise
                                    ConvergenceChecker = convergence_checker_FS();
                            end
                            NLSolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                    end
                    NLSolver.MaxIter = obj.SimulatorSettings.MaxIterations;
                    ConvergenceChecker.Tol = obj.SimulatorSettings.Tolerance;
                    switch (obj.SimulatorSettings.Formulation)
                        case('Immiscible')
                            ConvergenceChecker.NormCalculator = norm_calculator_immiscible();
                        otherwise
                            ConvergenceChecker.NormCalculator = norm_calculator_comp();
                    end
                    NLSolver.AddConvergenceChecker(ConvergenceChecker);
                    % Build FIM Coupling strategy
                    if obj.SimulatorSettings.LTS
                        Coupling = FIM_Strategy_LTS('FIM', NLSolver);
                    else
                        Coupling = FIM_Strategy('FIM', NLSolver);
                    end
                
                case('Sequential')
                    % Sequential coupling
                    Coupling = Sequential_Strategy('Sequential');
                    Coupling.MaxIter = obj.SimulatorSettings.MaxIterations;
                    % pressuresolver = incompressible_pressure_solver();
                    pressuresolver = NL_Solver();
                    pressuresolver.MaxIter = 15;
                    ConvergenceChecker = convergence_checker_pressure();
                    ConvergenceChecker.Tol = 1e-6;
                    pressuresolver.AddConvergenceChecker(ConvergenceChecker);
                    pressuresolver.SystemBuilder = pressure_system_builder();
                    switch obj.SimulatorSettings.DiscretizationModel
                        case ('ADM')
                            pressuresolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            pressuresolver.LinearSolver.OperatorsAssembler = operators_assembler_seq(1, 1);
                            pressuresolver.ConvergenceChecker.adm = 1;
                            pressuresolver.ConvergenceChecker.OperatorsAssembler = operators_assembler_seq(1,1);
                        case('MMs')
                            pressuresolver.LinearSolver = linear_solver_MMs(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            pressuresolver.LinearSolver.MSFE = MMsSettings.MSFE;
                            if MMsSettings.CorrectionFunctions
                                pressuresolver.LinearSolver.CorrectionFunctions = true;
                            end
                        otherwise
                            pressuresolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                    end
                    if obj.incompressible
                        pressuresolver.ConvergenceChecker.Incompressible = 1;
                    end
                    Coupling.AddPressureSolver(pressuresolver);
                    switch (obj.SimulatorSettings.TransportSolver.Type)
                        case('IMPSAT') 
                            transportsolver = NL_Solver();
                            transportsolver.MaxIter = obj.SimulatorSettings.TransportSolver.MaxIter;
                            ConvergenceChecker = convergence_checker_transport();
                            ConvergenceChecker.Tol = obj.SimulatorSettings.TransportSolver.Tol;
                            transportsolver.AddConvergenceChecker(ConvergenceChecker);
                            transportsolver.SystemBuilder = transport_system_builder();
                            Coupling.ConvergenceChecker = convergence_checker_outer();
                            Coupling.ConvergenceChecker.Tol = obj.SimulatorSettings.Tolerance;
                            switch obj.SimulatorSettings.DiscretizationModel
                                case('ADM')
                                    transportsolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                                    transportsolver.LinearSolver.OperatorsAssembler = operators_assembler_seq(2, obj.NofEq);
                                    transportsolver.ConvergenceChecker.adm = 1;
                                    transportsolver.ConvergenceChecker.OperatorsAssembler = operators_assembler_seq(2, obj.NofEq);
                                otherwise
                                    transportsolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            end
                        otherwise
                            transportsolver = explicit_transport_solver();
                            transportsolver.SystemBuilder = explicit_transport_system_builder();
                            Coupling.ConvergenceChecker = convergence_checker_impes();
                    end
                    Coupling.AddTransportSolver(transportsolver);
                case('SinglePhase')
                    % Single phase coupling
                    Coupling = SinglePhase_Strategy('SinglePhase');
                    pressuresolver = NL_Solver();
                    pressuresolver.MaxIter = 15;
                    ConvergenceChecker = convergence_checker_pressure();
                    ConvergenceChecker.Tol = 1e-6;
                    pressuresolver.AddConvergenceChecker(ConvergenceChecker);
                    pressuresolver.SystemBuilder = pressure_system_builder();
                    switch obj.SimulatorSettings.DiscretizationModel
                        case('ADM')
                            pressuresolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            pressuresolver.LinearSolver.OperatorsAssembler = operators_assembler_fim(obj.NofEq);
                            pressuresolver.ConvergenceChecker.adm = 1;
                            pressuresolver.ConvergenceChecker.OperatorsAssembler = operators_assembler_fim(obj.NofEq);
                        case ('MMs')
                            pressuresolver.LinearSolver = linear_solver_MMs(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            pressuresolver.LinearSolver.OperatorsAssembler = operators_assembler_seq(1, 1);
                            pressuresolver.LinearSolver.MSFE = obj.SimulatorSettings.MMsSettings.MSFE;
                            if obj.SimulatorSettings.MMsSettings.CorrectionFunctions
                                pressuresolver.LinearSolver.CorrectionFunctions = true;
                            end
                        otherwise
                            pressuresolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                    end
                    if obj.incompressible
                        Coupling.Incompressible = 1;
                        pressuresolver.ConvergenceChecker.Incompressible = 1;
                    end
                    Coupling.AddPressureSolver(pressuresolver);
            end
            Coupling.TimeStepSelector = timestep_selector(obj.SimulatorSettings.cfl, obj.SimulatorSettings.MinMaxdt(1), obj.SimulatorSettings.MinMaxdt(2));
            TimeDriver.AddCouplingStrategy(Coupling);
            switch(obj.SimulatorSettings.StopCriterion)
                case('MAX TIME')
                    end_of_sim_eval = end_of_sim_evaluator(obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps);
                case('COMPONENT CUT')
                    end_of_sim_eval = end_of_sim_evaluator_gascut(obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps, 0.2);
                case('PV INJECTED')
                    end_of_sim_eval = end_of_sim_evaluator_PVInjected(obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps, 3);
            end
            TimeDriver.AddEndOfSimEvaluator(end_of_sim_eval);
    end
        function Summary = BuildSummary(obj, simulation)
            %%%%%%%%%%%%%%% BuildObjects for OUTPUT%%%%%%%%%
            switch(obj.SimulatorSettings.CouplingType)
                case('FIM')
                    CouplingStats = FIM_Stats(obj.SimulatorSettings.MaxNumTimeSteps);
                case('Sequential')
                    CouplingStats = Sequential_Stats(obj.SimulatorSettings.MaxNumTimeSteps);
                case('SinglePhase')
                    CouplingStats = SinglePhase_Stats(obj.SimulatorSettings.MaxNumTimeSteps);
            end
            wellsData = wells_data(obj.SimulatorSettings.MaxNumTimeSteps, simulation.FluidModel.NofPhases, simulation.FluidModel.NofComp, simulation.ProductionSystem.Wells);
            switch (obj.SimulatorSettings.DiscretizationModel)
                case('ADM')
                    Summary = Run_Summary_ADM(obj.SimulatorSettings.MaxNumTimeSteps, CouplingStats, wellsData, simulation.DiscretizationModel.maxLevel(1)); % Only reservoir for now
                otherwise
                    Summary = Run_Summary(obj.SimulatorSettings.MaxNumTimeSteps, CouplingStats, wellsData);
            end
        end
        function Writer = BuildWriter(obj, InputDirectory, simulation)
            % Build Plotter
            switch(obj.SimulatorSettings.plotting)
                case('Matlab')
                    if simulation.DiscretizationModel.ReservoirGrid.Nx == 1 || simulation.DiscretizationModel.ReservoirGrid.Ny == 1
                        plotter = Matlab_Plotter_1D();
                    else
                        plotter = Matlab_Plotter_2D();
                    end
                case('VTK')
                    plotter = VTK_Plotter(InputDirectory, obj.SimulationInput.ProblemName);
                otherwise
                    warning('WARNING: NO valid Plotter was selected. Results will not be plotted.');
                    plotter = no_Plotter();
            end
            
            switch(obj.SimulatorSettings.DiscretizationModel)
                case ('ADM')
                    Writer = output_writer_adm(InputDirectory, obj.SimulationInput.ProblemName,...
                        simulation.ProductionSystem.Wells.NofInj, simulation.ProductionSystem.Wells.NofProd, ...
                        simulation.Summary.CouplingStats.NTimers, simulation.Summary.CouplingStats.NStats,...
                        simulation.FluidModel.NofComp);
                otherwise
                    Writer = output_writer_FS(InputDirectory, obj.SimulationInput.ProblemName,...
                        simulation.ProductionSystem.Wells.NofInj, simulation.ProductionSystem.Wells.NofProd,...
                        simulation.Summary.CouplingStats.NTimers, simulation.Summary.CouplingStats.NStats,...
                        simulation.FluidModel.NofComp);   
            end
            Writer.AddPlotter(plotter);
        end
        function DefineProperties(obj, ProductionSystem, FluidModel, DiscretizationModel)
            switch(obj.SimulationInput.FracturesProperties.Fractured)
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