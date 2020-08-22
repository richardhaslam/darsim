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
        capillarymodel = 0;
        NofEq
    end
    methods
        function simulation = BuildSimulation(obj, FractureMatrix)
            simulation = Reservoir_Simulation();
            simulation.DiscretizationModel = obj.BuildDiscretization(FractureMatrix);
            simulation.ProductionSystem = obj.BuildProductionSystem(FractureMatrix,simulation.DiscretizationModel);
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
            VarValues = ones(N, size(obj.SimulationInput.InitData,2));
            for i=1:size(obj.SimulationInput.InitData,2)
                VarValues(:, i) = VarValues(:, i) .* obj.SimulationInput.InitData(:,i);
            end        
            
            switch(simulation.FluidModel.name)
                case('SinglePhase')
                    VarNames = {'P_1', 'S_1'};
                    VarValues(:, 2) = 1;
                    simulation.Initializer = initializer_singlephase(VarNames, VarValues);
                case('Immiscible')
                    VarNames = {'P_2', 'S_1', 'S_2'};
                    VarValues(:, 3) = 1 - VarValues(:, 2);  % 3rd column is S2 which is 1-S1
                    simulation.Initializer = initializer(VarNames, VarValues);
                case('Geothermal_SinglePhase')
                    VarNames = {'P_1', 'T', 'S_1'};
                    VarValues(:, 3) = 1;
                    VarValues(:, 2) = obj.SimulationInput.ReservoirProperties.Temperature;
                    simulation.Initializer = initializer_singlephase(VarNames, VarValues);
                case('Geothermal_MultiPhase')
                    error('Geothermal_MultiPhase is not yet implemented');
                otherwise
                    VarNames = {'P_2', 'z_1', 'z_2'};
                    simulation.Initializer = initializer_hydrostatic(VarNames, VarValues);
            end
        end
        function Discretization = BuildDiscretization(obj, FractureMatrix)
            %% 1. Create fine-scale grids
            % 1. Reservoir Grid
            switch obj.SimulationInput.ReservoirProperties.Discretization
                case('CartesianGrid')
                    ReservoirGrid = cartesian_grid(obj.SimulationInput.ReservoirProperties.Grid.N);
                case('CornerPointGrid')
                    ReservoirGrid = corner_point_grid(obj.SimulationInput.ReservoirProperties);
                otherwise
                    error('At this moment, only "CartesianGrid" and "CornerPointGrid" discretization models are supported in DARSim!\n');
            end
            if obj.SimulationInput.FracturesProperties.isFractured
                [FracturesGrid, CrossConnections] = obj.ScanFracturesData(FractureMatrix, ReservoirGrid);
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
                                switch obj.SimulatorSettings.Formulation
                                    case {'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                                        prolongationbuilder.BFUpdater = bf_updater_ms_geothermal();
                                    otherwise
                                        prolongationbuilder.BFUpdater = bf_updater_ms();
                                end
                            else % Fractured
                                switch obj.SimulatorSettings.Formulation
                                    case {'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                                        prolongationbuilder.BFUpdater = bf_updater_FAMS_geothermal();
                                    otherwise
                                        prolongationbuilder.BFUpdater = bf_updater_FAMS();
                                end
                                prolongationbuilder.BFUpdater.BFtype = ADMSettings.BFtype;
                            end
                            if strcmp(ADMSettings.PInterpolator, 'Homogeneous')
                                prolongationbuilder.BFUpdater.MaxContrast = 1;
                            else
                               prolongationbuilder.BFUpdater.MaxContrast = ADMSettings.BF_MaxContrast;
                            end
                            prolongationbuilder.BFUpdater.pEDFM_MaxContrast = ADMSettings.pEDFM_MaxContrast;
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
                                prolongationbuilder.key = strcat('S_',num2str(i-1));
                        end
                        operatorshandler.AddProlongationBuilder(prolongationbuilder, i);
                    end
                    if obj.SimulatorSettings.Formulation == "Geothermal_SinglePhase" || obj.SimulatorSettings.Formulation == "Geothermal_MultiPhase"
                        switch(ADMSettings.HInterpolator)
                            case('Constant')
                                prolongationbuilder = prolongation_builder_constant(ADMSettings.maxLevel(1));
                            case('MS')
                                prolongationbuilder = prolongation_builder_MSHyperbolic(ADMSettings.maxLevel(1));
                                prolongationbuilder.key = 'T';
                        end
                        i = length(operatorshandler.ProlongationBuilders);
                        operatorshandler.AddProlongationBuilder(prolongationbuilder, i+1);
                    end
                    
                    % b. Grid selection criterion (time\space based)
                    switch (ADMSettings.GridSelCriterion)
                        case('dfdx')
                            gridselector = adm_grid_selector_delta(ADMSettings.tol, ADMSettings.key);
                        case('dfdx2')
                            gridselector = adm_grid_selector_delta2(ADMSettings.tol, obj.SimulatorSettings.LTStol, ADMSettings.key);
                        case('dfdt')
                            gridselector = adm_grid_selector_time(ADMSettings.tol, ADMSettings.key, obj.SimulationInput.ReservoirProperties.Grid.N_ActiveCells, ADMSettings.maxLevel(1));
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
                        prolongationbuilder.BFUpdater.BFtype = MMsSettings.BFtype;
                    end
                    % Reduce contrast for BF computation to remove peaks
                    prolongationbuilder.BFUpdater.MaxContrast = 10^-2;
                    prolongationbuilder.BFUpdater.MaxContrast = MMsSettings.BF_MaxContrast;
                    prolongationbuilder.BFUpdater.pEDFM_MaxContrast = MMsSettings.pEDFM_MaxContrast;
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
                    switch obj.SimulationInput.ReservoirProperties.Discretization
                        case('CartesianGrid')
                            Discretization = Cartesian_Discretization_model();
                        case('CornerPointGrid')
                            Discretization = CornerPointGrid_Discretization_model();
                            Discretization.CornerPointGridData = obj.SimulationInput.ReservoirProperties.CornerPointGridData;
                        otherwise
                            error('At this moment, only "CartesianGrid" and "CornerPointGrid" discretization models are supported in DARSim!\n');
                    end
            end
            
            %% 3. Add Grids to the Discretization Model
            Discretization.AddReservoirGrid(ReservoirGrid);
            if obj.SimulationInput.FracturesProperties.isFractured
                Discretization.AddFracturesGrid(FracturesGrid);
                Discretization.AddCrossConnections(CrossConnections);
            end
        end
        function ProductionSystem = BuildProductionSystem (obj, FractureMatrix, DiscretizationModel)
            ProductionSystem = Production_System();
            %% RESERVOIR
            Lx = obj.SimulationInput.ReservoirProperties.size(1);       %Dimension in x-direction [m]
            Ly = obj.SimulationInput.ReservoirProperties.size(2);       %Dimension in y-direction [m]
            h  = obj.SimulationInput.ReservoirProperties.size(3);       %Reservoir thickness (z-direction) [m]
            Tres = obj.SimulationInput.ReservoirProperties.Temperature; %Res temperature [K]
            Reservoir = reservoir(Lx, Ly, h, Tres);
            cr = obj.SimulationInput.ReservoirProperties.Compressibility;
            nx = obj.SimulationInput.ReservoirProperties.Grid.N(1);
            ny = obj.SimulationInput.ReservoirProperties.Grid.N(2);
            nz = obj.SimulationInput.ReservoirProperties.Grid.N(3);
            N_ActiveCells = obj.SimulationInput.ReservoirProperties.Grid.N_ActiveCells;
            Cpr = obj.SimulationInput.ReservoirProperties.SpecificHeat;
            RockDensity = obj.SimulationInput.ReservoirProperties.RockDensity;
            % Permeability Data
            K = ones(N_ActiveCells, 3);
            switch obj.SimulationInput.ReservoirProperties.PermUnit
                case('m2')
                    PermMultiplier = 1;
                case('D')
                    PermMultiplier = 1e-12;
                case('mD')
                    PermMultiplier = 1e-15;
            end
            if ~isempty(obj.SimulationInput.ReservoirProperties.Perm) && length(obj.SimulationInput.ReservoirProperties.Perm) == N_ActiveCells
                K = obj.SimulationInput.ReservoirProperties.Perm;
            else
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
                        if (nx~=field(1,1)) || (ny~=field(2,1)) || (nz~=field(3,1))
                            warning('The grid cells mentioned in the permeability file #%d do not match with Reservoir grid cells.\n',i);
                        end
                        if (nx>field(1,1)) || (ny>field(2,1)) || (nz>field(3,1))
                            error('The grid cells mentioned in the permeability file #%d are less than that of Reservoir grid cells. Check the input file.\n',i);
                        end
                        field1 = reshape(field(4:end,1),[field(1,1) field(2,1) field(3,1)]);
                        % make it the size of the grid
                        K(:,i) = reshape(field1(1:nx, 1:ny, 1:nz)*1e-15, nx*ny*nz, 1);
                        if strcmp(obj.SimulationInput.ReservoirProperties.PermScale, 'Logarithmic')
                            % In case the data is in logarithmic scale
                            K(:,i) = reshape(10.^(field1(1:nx, 1:ny, 1:nz)) * PermMultiplier, nx*ny*nz, 1);
                        else
                            K(:,i) = reshape(field1(1:nx, 1:ny, 1:nz) * PermMultiplier, nx*ny*nz, 1);
                        end
                        
                        % Reduce the heterogeneity contrast
                        if obj.SimulationInput.ReservoirProperties.PermContrastReduction
                            K_Log10 = log10(K(:,i));
                            if strcmp(obj.SimulationInput.ReservoirProperties.PermContrastMean,'Default')
                                K_Mean = mean(K(:, i));
                            else
                                K_Mean = str2double(obj.SimulationInput.ReservoirProperties.PermContrastMean) * PermMultiplier;
                            end
                            maxContrast = obj.SimulationInput.ReservoirProperties.PermContrastOrder;
                            
                            ratio = (max(K_Log10) - min(K_Log10)) / maxContrast;
                            K_Log10 = (  ( K_Log10-log10(K_Mean) ) / ratio  )  +  log10(K_Mean);
                            K(:, i) = 10.^(K_Log10);
                        end
                        
                    else
                        % Homogeneous Permeability
                        value = obj.SimulationInput.ReservoirProperties.Perm(i);
                        if strcmp(obj.SimulationInput.ReservoirProperties.PermScale, 'Logarithmic')
                            K(:, i)= K(:,i) * 10^value * PermMultiplier;
                        else
                            K(:, i)= K(:,i) * value * PermMultiplier;
                        end
                    end
                end
            end
            
            %%%% Adding Porosity Data
            if obj.SimulationInput.ReservoirProperties.PorosityInclude
                fprintf('\n---> Reading porosity file ...');
                field = load(obj.SimulationInput.ReservoirProperties.PorosityFile);
                fprintf(' ---> Completed.\n');
                
                % reshape it to specified size
                if (nx~=field(1,1)) || (ny~=field(2,1)) || (nz~=field(3,1))
                    warning('The grid cells mentioned in the porosity file do not match with Reservoir grid cells.\n');
                end
                if (nx>field(1,1)) || (ny>field(2,1)) || (nz>field(3,1))
                    error('The grid cells mentioned in the porosity file are less than that of Reservoir grid cells. Check the input file.\n');
                end
                field1 = reshape(field(4:end,1),[field(1,1) field(2,1) field(3,1)]);
                % make it the size of the grid
                phi = reshape(field1(1:nx, 1:ny, 1:nz), nx*ny*nz, 1);
            else
                phi = obj.SimulationInput.ReservoirProperties.phi .* ones(N_ActiveCells, 1);
            end
            if any(phi>1) || any(phi<0)
                error('In the porosity data, there are values out of the physical range. Please check the Input file!')
            end
            
            %%%% Adding permeability and Porosity info to the reservoir
            Reservoir.AddPermeabilityPorosity(K, phi);

            % Adding heat conductivity info to the reservoir
            if contains(obj.SimulatorSettings.Formulation,'Geothermal')
                Reservoir.AddConductivity(obj.SimulationInput.ReservoirProperties.RockConductivity,obj.SimulationInput.FluidProperties.FluidConductivity);
            end
            
            Reservoir.Cr = cr;
            Reservoir.Cpr = Cpr;
            Reservoir.Rho = RockDensity;
            Reservoir.P0 = obj.SimulationInput.InitData(:,1); % Initial Pressure of the Reservoir
            
            switch obj.SimulatorSettings.DiscretizationModel
                case('ADM')
                    % This is for DLGR type ADM: it reads coarse permeabilities
                    if obj.SimulatorSettings.ADMSettings.DLGR
                        K_coarse = cell(obj.SimulatorSettings.ADMSettings.maxLevel + 1, 1);
                        K_coarse{1} = K;
                        for L= 2:obj.SimulatorSettings.ADMSettings.maxLevel + 1
                            for d=1:2
                                % load the file in a vector
                                field = load(obj.SimulationInput.ReservoirProperties.CoarsePermFile{L-1,d});
                                % reshape it to specified size
                                k = field(4:end)*1e-15; % for now the cparse permeabilities are in [mD] unit
                                K_coarse{L}(:, d) = k;
                            end
                            K_coarse{L}(:, 3) = k;
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
            switch obj.SimulationInput.ReservoirProperties.Discretization
                case('CartesianGrid')
                    GridVolume = dx*dy*dz * ones(N_ActiveCells,1);
                case('CornerPointGrid')
                    GridVolume = DiscretizationModel.CornerPointGridData.Cells.Volume;
                otherwise
            end
            
            %Injectors
            for i=1:Wells.NofInj
                switch(obj.SimulationInput.WellsInfo.Inj(i).Formula.type)
                    case ('DIRICHLET')
                        PI = dy*dz/(dx/2) .* ones(N_ActiveCells,1);
                    case ('PI')
                        PI = obj.SimulationInput.WellsInfo.Inj(i).Formula.value .* ones(N_ActiveCells,1);
                    case('WI')
                        PI = obj.SimulationInput.WellsInfo.Inj(i).Formula.value .* GridVolume;
                    case ('RADIUS')
                        radius = obj.SimulationInput.WellsInfo.Inj(i).Formula.value;
                        error('DARSim Error: Radius calculation of PI is not implemented for now.')
                end
                Coordinate = obj.SimulationInput.WellsInfo.Inj(i).Coordinate;
                switch (obj.SimulationInput.FluidProperties.FluidModel)
                    case {'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                        temperature = obj.SimulationInput.WellsInfo.Inj(i).Temperature;
                    otherwise
                        temperature = Tres;
                end
                switch (obj.SimulationInput.WellsInfo.Inj(i).Constraint.name)
                    case('PRESSURE')
                        pressure = obj.SimulationInput.WellsInfo.Inj(i).Constraint.value;
                        % temperature = obj.SimulationInput.WellsInfo.Inj(i).Temperature;
                        Injector = injector_pressure(PI, Coordinate, pressure, temperature, n_phases);
                    case('RATE')
                        rate = obj.SimulationInput.WellsInfo.Inj(i).Constraint.value;
                        p_init = obj.SimulationInput.InitData(:,1);
                        rate = rate * Reservoir.TotalPV / (3600 * 24); % convert pv/day to m^3/s
                        Injector = injector_rate(PI, Coordinate, rate, p_init, temperature, n_phases);
                end
                Wells.AddInjector(Injector);
            end
            
            %Producers
            for i=1:Wells.NofProd
                switch(obj.SimulationInput.WellsInfo.Prod(i).Formula.type)
                    case ('DIRICHLET')
                        PI = dy*dz/(dx/2) .* ones(N_ActiveCells,1);
                    case ('PI')
                        PI = obj.SimulationInput.WellsInfo.Prod(i).Formula.value .* ones(N_ActiveCells,1);
                    case('WI')
                        PI = obj.SimulationInput.WellsInfo.Prod(i).Formula.value .* GridVolume;
                    case ('RADIUS')
                        radius = obj.SimulationInput.WellsInfo.Prod(i).Formula.value;
                        error('DARSim2 Error: Radius calculation of PI is not implemented for now.')
                end
                Coordinate = obj.SimulationInput.WellsInfo.Prod(i).Coordinate;
                switch (obj.SimulationInput.WellsInfo.Prod(i).Constraint.name)
                    case('PRESSURE')
                        pressure = obj.SimulationInput.WellsInfo.Prod(i).Constraint.value;
                        Producer = producer_pressure(PI, Coordinate, pressure);
                    case('RATE')
                        rate = obj.SimulationInput.WellsInfo.Prod(i).Constraint.value;
                        rate = rate * Reservoir.TotalPV / (3600 * 24); % convert pv/day to m^3/s
                        Producer = producer_rate(PI, Coordinate, rate);
                end
                Wells.AddProducer(Producer);
            end
            ProductionSystem.AddWells(Wells);
            
            %% FRACTURE NETWORK
            if obj.SimulationInput.FracturesProperties.isFractured
                FracturesNetwork = fracture_system();
                FracturesNetwork.Active = 1;
                temp = strfind(FractureMatrix, 'NUM_FRACS');
                index = find(~cellfun('isempty', temp));
                temp = strsplit(FractureMatrix{index},' ');
                NumOfFrac = temp{end};
                FracturesNetwork.NumOfFrac = str2double( NumOfFrac );
                
                temp = strfind(FractureMatrix, 'PROPERTIES');
                frac_index = find(~cellfun('isempty', temp));
                
                FracturesNetwork.Fractures = fracture();
                for f = 1 : FracturesNetwork.NumOfFrac
                    frac_info_split = strsplit(FractureMatrix{frac_index(f)},' ');                      % Splitted data for each fracture
                    FracturesNetwork.Fractures(f).Length = str2double( frac_info_split{3} );            % Length of each fracture
                    FracturesNetwork.Fractures(f).Width = str2double( frac_info_split{4} );             % Width of each fracture
                    FracturesNetwork.Fractures(f).Thickness = str2double( frac_info_split{5} );         % Thickness of each fracture
                    FracturesNetwork.Fractures(f).Temp = Tres;                                          % Temperature of each fracture
                    FracturesNetwork.Fractures(f).Cpr = Cpr;
                    FracturesNetwork.Fractures(f).Rho = RockDensity;
                    
                    Porosity = str2double( frac_info_split{6} );                                        % Porosity  of each fracture
                    Permeability = str2double( frac_info_split{7} );                                    % Permeability of each fracture
                    Kx = ones(DiscretizationModel.FracturesGrid.N(f), 1)*Permeability;
                    Ky = Kx;
                    Kz = Kx;
                    K = [Kx, Ky, Kz];
                    FracturesNetwork.Fractures(f).AddPermeabilityPorosity(K, Porosity);                 % Adding porosity and permeability to the fracture
                    if contains(obj.SimulatorSettings.Formulation,'Geothermal')
                        FracturesNetwork.Fractures(f).AddConductivity(obj.SimulationInput.ReservoirProperties.RockConductivity,obj.SimulationInput.FluidProperties.FluidConductivity);
                    end
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
                    n_comp = obj.SimulationInput.FluidProperties.NofComponents;
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
                case{'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                    % build the geothermal fluid model
                    switch (obj.SimulationInput.FluidProperties.FluidModel)
                        case{'Geothermal_SinglePhase'}
                            FluidModel = geothermal_singlephase_fluid_model();
                        case{'Geothermal_MultiPhase'}
                            FluidModel = geothermal_multiphase_fluid_model();
                    end
                    Phase = geothermal_singlephase_phase();
                    FluidModel.AddPhase(Phase, 1);
                    obj.SimulatorSettings.CouplingType = 'FIM';
                    %Gets all densities [kg/m^3]
                    Phase.rho0 = obj.SimulationInput.FluidProperties.Density;
                    %Gets all viscosities [Pa sec]
                    Phase.mu0 = obj.SimulationInput.FluidProperties.mu;
                    % Compressibility
                    Phase.cf0 = obj.SimulationInput.FluidProperties.Comp;
                    % Conductivity
                    Phase.Kf = obj.SimulationInput.FluidProperties.FluidConductivity;
                    % Specific Heat
                    Phase.Cp = obj.SimulationInput.FluidProperties.SpecificHeat;
                    FluidModel.AddPhase(Phase, 1);
                    if Phase.cf0 == 0
                        obj.incompressible = 1;
                    end
                otherwise
                    error('FluidModel is not defined!');
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
                case('Table')
                    FluidModel.RelPermModel = relperm_model_table(obj.SimulationInput.FluidProperties.RelPerm.TableType,...
                                                                  obj.SimulationInput.FluidProperties.RelPerm.TableData);
                    obj.SimulationInput.FluidProperties.RelPerm.s_irr = FluidModel.RelPermModel.S_irr;
            end
            
            % Irriducible sat
            for i=1:FluidModel.NofPhases
                FluidModel.Phases(i).sr = obj.SimulationInput.FluidProperties.RelPerm.s_irr(i);
            end
            %% Capillary pressure model
            switch (obj.SimulationInput.FluidProperties.Capillarity.name)
                case('JLeverett')
                    obj.capillarymodel = 1;
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
                    if obj.SimulatorSettings.LTS
                        Formulation.MatrixAssembler = LTS_matrix_assembler();
                    else
                        Formulation.MatrixAssembler = matrix_assembler();
                    end
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
                case('Geothermal_SinglePhase')
                    Formulation = geothermal_singlephase_formulation();
                    Formulation.MatrixAssembler = matrix_assembler_geothermal();
                    obj.NofEq = 2;
                case('Geothermal_MultiPhase')
                    Formulation = geothermal_multiphase_formulation();
                    obj.NofEq = 2;
            end
            Formulation.NofPhases = obj.SimulationInput.FluidProperties.NofPhases;
        end
        function TimeDriver = BuildTimeDriver(obj)
            if obj.SimulatorSettings.LTSPlot == 1
                TimeDriver = LTSPrint_TimeLoop_Driver(obj.SimulatorSettings.Reports, obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps);
                % elseif obj.SimulatorSettings.ADT_SEQ == 1
                %   TimeDriver = LTS_TimeLoop_Driver(obj.SimulatorSettings.Reports, obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps);
            else
                TimeDriver = TimeLoop_Driver(obj.SimulatorSettings.Reports, obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps);
            end
            % Construct Coupling
            switch(obj.SimulatorSettings.CouplingType)
                case('FIM')
                    % FIM coupling
                    %%%%FIM settings
                    switch obj.SimulatorSettings.Formulation
                        case {'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                            NLSolver = NL_Solver_geothermal();
                        otherwise
                            NLSolver = NL_Solver();
                    end
                    NLSolver.SystemBuilder = fim_system_builder();
                    switch obj.SimulatorSettings.DiscretizationModel
                        case ('ADM')
                            % Build a different convergence cheker and a proper LS for ADM
                            switch obj.SimulatorSettings.Formulation
                                case {'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                                    ConvergenceChecker = convergence_checker_ADM_Geothermal();
                                otherwise
                                    ConvergenceChecker = convergence_checker_ADM();
                            end
                            
                            ConvergenceChecker.OperatorsAssembler = operators_assembler_fim(obj.NofEq);
                            NLSolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            NLSolver.LinearSolver.OperatorsAssembler = operators_assembler_fim(obj.NofEq);
                            if obj.SimulatorSettings.ADMSettings.DLGR
                                % it will change perm during ADM
                                % simulaiton to use upscaled ones
                                NLSolver.LinearSolver.DLGR = 1;
                            end
                        otherwise
                            switch (obj.SimulatorSettings.Formulation)
                                case('Molar')
                                    ConvergenceChecker = convergence_checker_finescale_molar();
                                case{'Geothermal_SinglePhase'}
                                    ConvergenceChecker = convergence_checker_finescale_geothermal_singlephase();
                                case{'Geothermal_MultiPhase'}
                                    ConvergenceChecker = convergence_checker_finescale_geothermal_multiphase();
                                otherwise
                                    ConvergenceChecker = convergence_checker_finescale();
                            end
                            NLSolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                    end
                    NLSolver.SystemBuilder.NumberOfEq = obj.NofEq;
                    NLSolver.MaxIter = obj.SimulatorSettings.MaxIterations;
                    ConvergenceChecker.NumberOfEq  = obj.NofEq;
                    if length(obj.SimulatorSettings.ResidualTolerances) < obj.NofEq
                        obj.SimulatorSettings.ResidualTolerances = obj.SimulatorSettings.ResidualTolerances(1) * ones(obj.NofEq, 1);
                    end
                    if length(obj.SimulatorSettings.SolutionTolerances) < obj.NofEq
                        obj.SimulatorSettings.SolutionTolerances = obj.SimulatorSettings.SolutionTolerances(1) * ones(obj.NofEq, 1);
                    end
                    ConvergenceChecker.ResidualTol = obj.SimulatorSettings.ResidualTolerances;
                    ConvergenceChecker.SolutionTol = obj.SimulatorSettings.SolutionTolerances;
                    switch (obj.SimulatorSettings.Formulation)
                        case('Immiscible')
                            ConvergenceChecker.NormCalculator = norm_calculator_immiscible();
                        case{'Geothermal_SinglePhase'}
                            ConvergenceChecker.NormCalculator = norm_calculator_geothermal_singlephase();
                        case{'Geothermal_MultiPhase'}
                            ConvergenceChecker.NormCalculator = norm_calculator_geothermal_multiphase();
                        otherwise
                            ConvergenceChecker.NormCalculator = norm_calculator_comp();
                    end
                    NLSolver.AddConvergenceChecker(ConvergenceChecker);
                    % Build FIM Coupling strategy
                    if obj.SimulatorSettings.LTS
                        LTSNLSolver = NL_Solver();
                        LTSNLSolver.SystemBuilder = LTS_fim_system_builder();
                        LTSNLSolver.SystemBuilder.NumberOfEq = obj.NofEq;
                        LTSNLSolver.MaxIter = obj.SimulatorSettings.MaxIterations;
                        LTSNLSolver.LinearSolver = LTS_linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                        LTSNLSolver.AddConvergenceChecker(ConvergenceChecker);
                        Coupling = LTS_FIM_Strategy('FIM', NLSolver, LTSNLSolver, obj.SimulatorSettings.LTStol);
                    else
                        Coupling = FIM_Strategy('FIM', NLSolver);
                    end
                    %% Sequential coupling
                case('Sequential')
                    % Sequential coupling
                    if obj.SimulatorSettings.LTS == 1 && strcmp(obj.SimulatorSettings.DiscretizationModel,'ADM') == 1
                        Coupling =  LTS_ADM_Adaptive_Sequential_Strategy('Sequential', obj.SimulatorSettings.LTStol);
                        
                        LTStransportsolver = LTS_NL_Solver();
                        LTStransportsolver.LinearSolver = LTS_linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                        
                        LTStransportsolver.MaxIter = obj.SimulatorSettings.TransportSolver.MaxIter;
                        ConvergenceChecker = convergence_checker_transport();
                        %ConvergenceChecker = LTS_convergence_checker_transport_ADM();
                        ConvergenceChecker.ResidualTol = obj.SimulatorSettings.TransportSolver.Tol;
                        ConvergenceChecker.SolutionTol = 1e-3;
                        LTStransportsolver.AddConvergenceChecker(ConvergenceChecker);
                        LTStransportsolver.ConvergenceChecker.adm = 1;
                        LTStransportsolver.ConvergenceChecker.lts = 1;
                        if obj.capillarymodel == 1
                            LTStransportsolver.ConvergenceChecker.capillary = 1;
                        end
                        
                        LTStransportsolver.SystemBuilder = LTStransport_system_builder();
                        LTStransportsolver.LinearSolver.OperatorsAssembler = LTS_operators_assembler_seq(2, obj.NofEq);
                        Coupling.AddLTSTransportSolver(LTStransportsolver);
                        
                        % OUTER ITERATION (for LTS and not)
                        Coupling.ConvergenceChecker = convergence_checker_outer();
                        Coupling.ConvergenceChecker.ResidualTol = obj.SimulatorSettings.TransportSolver.Tol;
                        Coupling.MaxIter = obj.SimulatorSettings.MaxIterations;
                        % pressuresolver = incompressible_pressure_solver();
                        pressuresolver = NL_Solver();
                        pressuresolver.MaxIter = 15;
                        ConvergenceChecker = convergence_checker_pressure();
                        ConvergenceChecker.ResidualTol = 1e-6;
                        pressuresolver.AddConvergenceChecker(ConvergenceChecker);
                        pressuresolver.SystemBuilder = pressure_system_builder();
                        
                        
                        % at the moment we use a fine pressure solver
                        %                     switch obj.SimulatorSettings.DiscretizationModel
                        %                         case ('ADM')
                        %                             pressuresolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                        %                             pressuresolver.LinearSolver.OperatorsAssembler = operators_assembler_seq(1, 1);
                        %                             pressuresolver.ConvergenceChecker.adm = 1;
                        %                             pressuresolver.ConvergenceChecker.OperatorsAssembler = operators_assembler_seq(1,1);
                        %                         case('MMs')
                        %                             pressuresolver.LinearSolver = linear_solver_MMs(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                        %                             MMsSettings = obj.SimulatorSettings.MMsSettings;
                        %                             pressuresolver.LinearSolver.MSFE = MMsSettings.MSFE;
                        %                             if MMsSettings.CorrectionFunctions
                        %                                 pressuresolver.LinearSolver.CorrectionFunctions = true;
                        %                             end
                        %                         otherwise
                        
                        pressuresolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                        
                        if obj.incompressible
                            pressuresolver.ConvergenceChecker.Incompressible = 1;
                        end
                        Coupling.AddPressureSolver(pressuresolver);
                        
                        transportsolver = NL_Solver();
                        transportsolver.MaxIter = obj.SimulatorSettings.TransportSolver.MaxIter;
                        ConvergenceChecker = convergence_checker_transport();
                        ConvergenceChecker.ResidualTol = obj.SimulatorSettings.TransportSolver.Tol;
                        ConvergenceChecker.SolutionTol = 0.2;
                        
                        
                        transportsolver.AddConvergenceChecker(ConvergenceChecker);
                        transportsolver.SystemBuilder = transport_system_builder();
                        Coupling.ConvergenceChecker = convergence_checker_outer();
                        Coupling.ConvergenceChecker.SolutionTol = obj.SimulatorSettings.Tolerance;
                        
                        transportsolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-2, 500);
                        transportsolver.LinearSolver.OperatorsAssembler = LTS_operators_assembler_seq(2, obj.NofEq);
                        transportsolver.ConvergenceChecker.adm = 1;
                        if obj.capillarymodel == 1
                            transportsolver.ConvergenceChecker.capillary = 1;
                        end
                        transportsolver.ConvergenceChecker.OperatorsAssembler = LTS_operators_assembler_seq(2, obj.NofEq);
                        
                        Coupling.AddTransportSolver(transportsolver);
                    else
                        if obj.SimulatorSettings.LTS == 1
                            switch (obj.SimulatorSettings.LTSCriterion)
                                case('Fixed')
                                    Coupling = LTS_Sequential_Strategy('Sequential', obj.SimulatorSettings.LTStol);
                                case('Adaptive')
                                    Coupling = LTS_Adaptive_Sequential_Strategy('Sequential', obj.SimulatorSettings.LTStol);
                            end
                            LTStransportsolver = LTS_NL_Solver();
                            LTStransportsolver.MaxIter = obj.SimulatorSettings.TransportSolver.MaxIter;
                            ConvergenceChecker = convergence_checker_transport();
                            ConvergenceChecker.ResidualTol = obj.SimulatorSettings.TransportSolver.Tol;
                            ConvergenceChecker.SolutionTol = 2e-2;
                            LTStransportsolver.AddConvergenceChecker(ConvergenceChecker);
                            if obj.capillarymodel == 1
                                LTStransportsolver.ConvergenceChecker.capillary = 1;
                            end
                            LTStransportsolver.SystemBuilder = LTStransport_system_builder();
                            LTStransportsolver.LinearSolver = LTS_linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            Coupling.AddLTSTransportSolver(LTStransportsolver);
                        else
                            Coupling = Sequential_Strategy('Sequential');
                        end
                        % OUTER ITERATION (for LTS and not)
                        Coupling.ConvergenceChecker = convergence_checker_outer();
                        Coupling.ConvergenceChecker.ResidualTol = obj.SimulatorSettings.TransportSolver.Tol;
                        Coupling.MaxIter = obj.SimulatorSettings.MaxIterations;
                        
                        pressuresolver = NL_Solver();
                        pressuresolver.MaxIter = 15;
                        ConvergenceChecker = convergence_checker_pressure();
                        ConvergenceChecker.ResidualTol = 1e-6;
                        pressuresolver.AddConvergenceChecker(ConvergenceChecker);
                        pressuresolver.SystemBuilder = pressure_system_builder();
                        switch obj.SimulatorSettings.DiscretizationModel
                            case ('ADM')
                                pressuresolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                                
                                %                                 pressuresolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                                %                                 pressuresolver.LinearSolver.OperatorsAssembler = operators_assembler_seq(1, 1);
                                %                                 pressuresolver.ConvergenceChecker.adm = 1;
                                %                                 pressuresolver.ConvergenceChecker.OperatorsAssembler = operators_assembler_seq(1,1);
                            case('MMs')
                                pressuresolver.LinearSolver = linear_solver_MMs(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                                MMsSettings = obj.SimulatorSettings.MMsSettings;
                                pressuresolver.LinearSolver.MSFE = MMsSettings.MSFE;
                                if MMsSettings.CorrectionFunctions
                                    pressuresolver.LinearSolver.CorrectionFunctions = true;
                                end
                            otherwise
                                pressuresolver.LinearSolver = linear_solver(obj.SimulatorSettings.LinearSolver, 1e-6, 1);
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
                                ConvergenceChecker.ResidualTol = obj.SimulatorSettings.TransportSolver.Tol;
                                ConvergenceChecker.SolutionTol = 1e-3;
                                transportsolver.AddConvergenceChecker(ConvergenceChecker);
                                if obj.capillarymodel == 1
                                    transportsolver.ConvergenceChecker.capillary = 1;
                                end
                                transportsolver.SystemBuilder = transport_system_builder();
                                Coupling.ConvergenceChecker = convergence_checker_outer();
                                Coupling.ConvergenceChecker.SolutionTol = obj.SimulatorSettings.Tolerance;
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
                    end
                case('SinglePhase')
                    % Single phase coupling
                    Coupling = SinglePhase_Strategy('SinglePhase');
                    pressuresolver = NL_Solver();
                    pressuresolver.MaxIter = 15;
                    ConvergenceChecker = convergence_checker_pressure();
                    ConvergenceChecker.ResidualTol = 1e-6;
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
            if obj.SimulatorSettings.LTS == 1 || obj.SimulatorSettings.FixedStep == 1
                Coupling.TimeStepSelector = fixed_timestep_selector(obj.SimulatorSettings.cfl, obj.SimulatorSettings.MinMaxdt(1), obj.SimulatorSettings.MinMaxdt(2));
            else
                Coupling.TimeStepSelector = timestep_selector(obj.SimulatorSettings.cfl, obj.SimulatorSettings.MinMaxdt(1), obj.SimulatorSettings.MinMaxdt(2));
            end
            Coupling.TimeStepSelector.TotalTime = obj.SimulationInput.TotalTime;
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
                    if obj.SimulatorSettings.LTS == 1
                        CouplingStats = LTS_FIM_Stats(obj.SimulatorSettings.MaxNumTimeSteps, 'LTS_FIM');
                    else
                        CouplingStats = FIM_Stats(obj.SimulatorSettings.MaxNumTimeSteps, 'FIM');
                    end
                    
                case('Sequential')
                    if obj.SimulatorSettings.LTS == 1
                        CouplingStats = LTS_Sequential_Stats(obj.SimulatorSettings.MaxNumTimeSteps, 'LTS_Sequential');
                    else
                        CouplingStats = Sequential_Stats(obj.SimulatorSettings.MaxNumTimeSteps, 'Sequential');
                    end
                case('SinglePhase')
                    CouplingStats = SinglePhase_Stats(obj.SimulatorSettings.MaxNumTimeSteps, 'SinglePhase');
            end
            wellsData = wells_data(obj.SimulatorSettings.MaxNumTimeSteps, simulation.FluidModel.NofPhases, simulation.FluidModel.NofComponents, simulation.ProductionSystem.Wells);
            switch (obj.SimulatorSettings.DiscretizationModel)
                case('ADM')
                    Summary = Run_Summary_ADM(obj.SimulatorSettings.MaxNumTimeSteps, CouplingStats, wellsData, simulation.DiscretizationModel.maxLevel(1)); % Only reservoir for now
                otherwise
                    Summary = Run_Summary(obj.SimulatorSettings.MaxNumTimeSteps, CouplingStats, wellsData);
            end
        end
        function Writer = BuildWriter(obj, InputDirectory, simulation)
            % Build Plotter
            switch(obj.SimulatorSettings.plotting.Software)
                case('Matlab')
                    if simulation.DiscretizationModel.ReservoirGrid.Nx == 1 || simulation.DiscretizationModel.ReservoirGrid.Ny == 1
                        plotter = Matlab_Plotter_1D();
                    else
                        plotter = Matlab_Plotter_2D();
                    end
                case({'ParaView','VTK'})
                    switch obj.SimulationInput.ReservoirProperties.Discretization
                        case('CornerPointGrid')
                            plotter = CornerPointGrid_VTK_Plotter(InputDirectory, obj.SimulationInput.ProblemName,obj.SimulatorSettings.NumOfPreviousReports);
                            plotter.PlotInterfaces = SimulatorSettings.plotting.PlotInterfaces;
                        otherwise
                            plotter = VTK_Plotter(InputDirectory, obj.SimulationInput.ProblemName,obj.SimulatorSettings.NumOfPreviousReports);
                    end
                    switch obj.SimulatorSettings.plotting.Format
                        case{'BINARY'}
                            plotter.isBinary = 1;
                        case{'ASCII'}    
                            plotter.isBinary = 0;
                        otherwise
                            warning('WARNING: NO valid output file format ("BINARY" or "ASCII") for Plotting was selected. Binary format is set by default.\n');
                            plotter.isBinary = 1;
                    end

                otherwise
                    warning('WARNING: NO valid Plotter was selected. Results will not be plotted.');
                    plotter = no_Plotter();
            end
            
            switch(obj.SimulatorSettings.DiscretizationModel)
                case ('ADM')
                    Writer = output_writer_adm(InputDirectory, obj.SimulationInput.ProblemName,...
                        simulation.ProductionSystem.Wells.NofInj,   simulation.ProductionSystem.Wells.NofProd, ...
                        simulation.Summary.CouplingStats.NTimers,   simulation.Summary.CouplingStats.NStats,...
                        obj.SimulatorSettings.NumOfPreviousReports, simulation.FluidModel.NofComp);
                otherwise
                    Writer = output_writer_FS(InputDirectory, obj.SimulationInput.ProblemName,...
                        simulation.ProductionSystem.Wells.NofInj,   simulation.ProductionSystem.Wells.NofProd,...
                        simulation.Summary.CouplingStats.NTimers,   simulation.Summary.CouplingStats.NStats,...
                        obj.SimulatorSettings.NumOfPreviousReports, simulation.FluidModel.NofComponents);
            end
            Writer.AddPlotter(plotter);
        end
        function DefineProperties(obj, ProductionSystem, FluidModel, DiscretizationModel)
            switch(obj.SimulationInput.FracturesProperties.isFractured)
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
        function [FracturesGrid, CrossConnections] = ScanFracturesData(obj, FractureMatrix, ReservoirGrid)
            fprintf('This simulation uses input from "%s" fracture generator.\n', obj.SimulationInput.FracturesProperties.FractureModelType );
            
            if ~strcmp( obj.SimulationInput.FracturesProperties.DiscretizationType, obj.SimulationInput.ReservoirProperties.Discretization )
                error('The discretization type used in "Fracture_Output" file (%s) does not match the one in the simulation input file (%s)!', ...
                        obj.SimulationInput.FracturesProperties.DiscretizationType, obj.SimulationInput.ReservoirProperties.Discretization );
            end
            
            NumOfFrac = obj.SimulationInput.FracturesProperties.NumOfFrac;
            fprintf('Extracting data from %02d fractures ...\n', NumOfFrac);
            FracturesGrid = fractures_grid(NumOfFrac);
            
            if ( obj.SimulationInput.ReservoirProperties.size(1) ~= obj.SimulationInput.FracturesProperties.Reservoir.Lx ) || ...
               ( obj.SimulationInput.ReservoirProperties.size(2) ~= obj.SimulationInput.FracturesProperties.Reservoir.Ly ) || ...
               ( obj.SimulationInput.ReservoirProperties.size(3) ~= obj.SimulationInput.FracturesProperties.Reservoir.Lz )
               warning('The dimension of reservoir in the "Fracture_Output" file does not match the simulation input file. This may cause inconsistency!');
            end
            
            if ( ReservoirGrid.Nx ~= obj.SimulationInput.FracturesProperties.Reservoir.Nx ) || ...
               ( ReservoirGrid.Ny ~= obj.SimulationInput.FracturesProperties.Reservoir.Ny ) || ...
               ( ReservoirGrid.Nz ~= obj.SimulationInput.FracturesProperties.Reservoir.Nz )
                error('The number of grid cells for reservoir in the "Fracture_Output" file does not match the simulation input file!');
            end

            temp = strfind(FractureMatrix, 'PROPERTIES');
            frac_index = find(~cellfun('isempty', temp));
            
            Nf = zeros(NumOfFrac,1);
            for f = 1 : NumOfFrac
                % Creat cartesian grid in each fracture
                frac_info_split = strsplit( FractureMatrix{frac_index(f)} , ' ' );
                grid_temp = strsplit(frac_info_split{8}, 'x');
                nx = str2double(grid_temp{1});
                ny = str2double(grid_temp{2});
                nz = 1;
                FractureGrid = cartesian_grid([nx;ny;nz]);
                
                % Add fracture grid coordinates (it's for plotting purposes)
                temp = strfind(FractureMatrix, 'GRID_COORDS_X');  frac_grid_coords_x = find(~cellfun('isempty', temp));  
                temp = strfind(FractureMatrix, 'GRID_COORDS_Y');  frac_grid_coords_y = find(~cellfun('isempty', temp));
                temp = strfind(FractureMatrix, 'GRID_COORDS_Z');  frac_grid_coords_z = find(~cellfun('isempty', temp));
                
                frac_grid_coords_x_split = str2double( strsplit(FractureMatrix{frac_grid_coords_x(f)},' ') );  frac_grid_coords_x_split(1) = [];
                frac_grid_coords_y_split = str2double( strsplit(FractureMatrix{frac_grid_coords_y(f)},' ') );  frac_grid_coords_y_split(1) = [];
                frac_grid_coords_z_split = str2double( strsplit(FractureMatrix{frac_grid_coords_z(f)},' ') );  frac_grid_coords_z_split(1) = [];
                FractureGrid.GridCoords = [ frac_grid_coords_x_split' , frac_grid_coords_y_split' , frac_grid_coords_z_split' ];
                FracturesGrid.AddGrid(FractureGrid, f);
                Nf(f) = nx*ny*nz;
            end
            
            % Reading the non-neighboring connectivities of frac-frac and frac-matrix
            NumOfAllFracGrids = obj.SimulationInput.FracturesProperties.NumOfAllFracGrids;
            CrossConnections = cross_connections();
            
            temp = strfind(FractureMatrix, '/');
            endOfSection = find(~cellfun('isempty', temp));
            
            temp = strfind(FractureMatrix, 'FRACCELL');
            fracCellIndeces = find(~cellfun('isempty', temp));
            fracCellIndeces = [fracCellIndeces; endOfSection(end)];
            
            temp = strfind(FractureMatrix, 'ROCK_CONN_EDFM');
            rockConnIndeces_EDFM = find(~cellfun('isempty', temp));
            
            temp = strfind(FractureMatrix, 'ROCK_CONN_pEDFM');
            rockConnIndeces_pEDFM = find(~cellfun('isempty', temp));
            
            temp = strfind(FractureMatrix, 'FRAC_CONN_EDFM');
            fracConnIndeces_EDFM = find(~cellfun('isempty', temp));
            
            temp = strfind(FractureMatrix, 'FRAC_CONN_pEDFM');
            fracConnIndeces_pEDFM = find(~cellfun('isempty', temp));
            
            % Number of phases (useful to define size of some objects)
            n_phases = obj.SimulationInput.FluidProperties.NofPhases;
            fprintf('---> Fracture ');
            
            Nm = ReservoirGrid.N;
            for f = 1 : NumOfFrac
                if (f>1),  fprintf(repmat('\b', 1, 9+27));  end
                fprintf('%04d/%04d',f,NumOfFrac);
                % looping over all global fracture cells
                fprintf(' ---> Grid cell ');
                for If = 1:Nf(f)
                    if (If>1),  fprintf(repmat('\b', 1, 11));  end
                    fprintf('%05d/%05d',If,Nf(f));
                    IfGlobal = sum(Nf(1:f-1))+If;
                    fracCellInfo_split = strsplit( FractureMatrix{fracCellIndeces(IfGlobal)} , {' ','	'} );
                    fracCellInd = str2double(fracCellInfo_split{2})+1;
                    if If~= fracCellInd,  error('The Indexing for fracture %02d cell %05d is not correct!',f,If);  end
                    NumOfRockConn_EDFM  = str2double(fracCellInfo_split{3});
                    NumOfRockConn_pEDFM = str2double(fracCellInfo_split{4});
                    NumOfFracConn_EDFM  = str2double(fracCellInfo_split{5});
                    NumOfFracConn_pEDFM = str2double(fracCellInfo_split{6});
                    
                    Counter=0;
                    
                    % Adding fracture-matrix EDFM connectivities
                    CrossConnections(IfGlobal,1).Cells = [];
                    rockConnInd_EDFM = rockConnIndeces_EDFM( rockConnIndeces_EDFM > fracCellIndeces(IfGlobal) );
                    rockConnInd_EDFM = rockConnInd_EDFM( rockConnInd_EDFM < fracCellIndeces(IfGlobal+1) );
                    if length(rockConnInd_EDFM) ~= NumOfRockConn_EDFM
                        error('The number of rock EDFM connectivities for fracture %02d cell %05d is not correct!',f,If);
                    end
                    for m = 1 : NumOfRockConn_EDFM
                        rockConnInfo = strsplit( FractureMatrix{rockConnInd_EDFM(m)} ,  {' ','	'} );
                        Im = str2double( rockConnInfo{2} ) + 1;
                        CI = str2double( rockConnInfo{3} );
                        Counter = Counter + 1;
                        CrossConnections( IfGlobal , 1 ).Cells( Counter ,1 ) = Im;
                        CrossConnections( IfGlobal , 1 ).CI   ( Counter ,1 ) = CI;
                    end
                    
                    % Adding fracture-matrix pEDFM connectivities
                    if strcmp(obj.SimulationInput.FracturesProperties.FractureModelType,'pEDFM')
                        rockConnInd_pEDFM = rockConnIndeces_pEDFM( rockConnIndeces_pEDFM > fracCellIndeces(IfGlobal) );
                        rockConnInd_pEDFM = rockConnInd_pEDFM( rockConnInd_pEDFM < fracCellIndeces(IfGlobal+1) );
                        if length(rockConnInd_pEDFM) ~= NumOfRockConn_pEDFM
                            error('The number of rock pEDFM connectivities for fracture %02d cell %05d is not correct!',f,If);
                        end
                        for m = 1 : NumOfRockConn_pEDFM
                            rockConnInfo = strsplit( FractureMatrix{rockConnInd_pEDFM(m)} ,  {' ','	'} );
                            Im = str2double( rockConnInfo{2} ) + 1;
                            CI = str2double( rockConnInfo{3} );
                            if ismember( Im, CrossConnections( IfGlobal , 1 ).Cells )
                                ind = find(CrossConnections( IfGlobal , 1 ).Cells == Im);
                                CrossConnections( IfGlobal , 1 ).CI( ind ,1 ) = CrossConnections( IfGlobal , 1 ).CI( ind ,1 ) + CI;
                            else
                                Counter = Counter + 1;
                                CrossConnections( IfGlobal , 1 ).Cells( Counter ,1 ) = Im;
                                CrossConnections( IfGlobal , 1 ).CI   ( Counter ,1 ) = CI;
                            end
                        end
                    end
                    
                    % Adding fracture-fracture EDFM connectivities
                    fracConnInd_EDFM = fracConnIndeces_EDFM( fracConnIndeces_EDFM > fracCellIndeces(IfGlobal) );
                    fracConnInd_EDFM = fracConnInd_EDFM( fracConnInd_EDFM < fracCellIndeces(IfGlobal+1) );
                    if length(fracConnInd_EDFM) ~= NumOfFracConn_EDFM
                        error('The number of frac EDFM connectivities for fracture %02d cell %05d is not correct!',f,If);
                    end
                    for m = 1 : NumOfFracConn_EDFM
                        fracConnInfo = strsplit( FractureMatrix{fracConnInd_EDFM(m)} ,  {' ','	'} );
                        g = str2double( fracConnInfo{2} ) + 1;
                        ig = str2double( fracConnInfo{3} ) + 1;
                        Ig = sum( Nf(1:g-1) )+ ig;
                        CI = str2double( fracConnInfo{4} );
                        if Ig > IfGlobal
                            Counter = Counter + 1;
                            CrossConnections( IfGlobal , 1 ).Cells( Counter ,1 ) = Nm+Ig;
                            CrossConnections( IfGlobal , 1 ).CI   ( Counter ,1 ) = CI;
                        end
                    end
                    
                    % Adding fracture-fracture pEDFM connectivities
                    if strcmp(obj.SimulationInput.FracturesProperties.FractureModelType,'pEDFM')
                        fracConnInd_pEDFM = fracConnIndeces_pEDFM( fracConnIndeces_pEDFM > fracCellIndeces(IfGlobal) );
                        fracConnInd_pEDFM = fracConnInd_pEDFM( fracConnInd_pEDFM < fracCellIndeces(IfGlobal+1) );
                        for m = 1 : NumOfFracConn_pEDFM
                            fracConnInfo = strsplit( FractureMatrix{fracConnInd_pEDFM(m)} ,  {' ','	'} );
                            g = str2double( fracConnInfo{2} ) + 1;
                            ig = str2double( fracConnInfo{3} ) + 1;
                            Ig = sum( Nf(1:g-1) )+ ig;
                            CI = str2double( fracConnInfo{4} );
                            if Ig > IfGlobal
                                if ismember( Nm+Ig, CrossConnections( IfGlobal , 1 ).Cells )
                                    ind = find(CrossConnections( IfGlobal , 1 ).Cells == Nm+Ig);
                                    CrossConnections( IfGlobal , 1 ).CI( ind ,1 ) = CrossConnections( IfGlobal , 1 ).CI( ind ,1 ) + CI;
                                else
                                    Counter = Counter + 1;
                                    CrossConnections( IfGlobal , 1 ).Cells( Counter ,1 ) = Nm+Ig;
                                    CrossConnections( IfGlobal , 1 ).CI   ( Counter ,1 ) = CI;
                                end
                            end
                        end
                    end
                    
                    % Initializing CrossConnections
                    [~, sort_ind] = sort( CrossConnections(IfGlobal,1).Cells );
                    CrossConnections(IfGlobal,1).Cells  = CrossConnections(IfGlobal,1).Cells(sort_ind);
                    CrossConnections(IfGlobal,1).CI     = CrossConnections(IfGlobal,1).CI(sort_ind);
                    CrossConnections(IfGlobal,1).UpWind = zeros(length(CrossConnections(IfGlobal,1).Cells), n_phases);
                    CrossConnections(IfGlobal,1).U_Geo  = zeros(length(CrossConnections(IfGlobal,1).Cells), n_phases);
                    CrossConnections(IfGlobal,1).T_Geo  = zeros(  size(CrossConnections(IfGlobal,1).CI)  );
                    CrossConnections(IfGlobal,1).T_Geo_Cond  = zeros(  size(CrossConnections(IfGlobal,1).CI)  );
                end
            end
            
            % Reading the pEDFM alpha corrections for reservoir transmissiblities
            if strcmp(obj.SimulationInput.FracturesProperties.FractureModelType,'pEDFM')
                
                % 1. Reservoir
                switch obj.SimulationInput.FracturesProperties.DiscretizationType
                    case('CartesianGrid')
                        pEDFM_alpha_Tx = zeros(ReservoirGrid.Nx+1,ReservoirGrid.Ny,ReservoirGrid.Nz);
                        pEDFM_alpha_Ty = zeros(ReservoirGrid.Nx,ReservoirGrid.Ny+1,ReservoirGrid.Nz);
                        pEDFM_alpha_Tz = zeros(ReservoirGrid.Nx,ReservoirGrid.Ny,ReservoirGrid.Nz+1);
                        % pEDFM_alpha_Tx
                        temp = strfind(FractureMatrix, 'ROCK_ALPHA_TX');
                        pEDFM_alpha_Tx_index = find(~cellfun('isempty', temp));
                        for t = 1 : length(pEDFM_alpha_Tx_index)
                            pEDFM_alpha_Tx_Split = strsplit(FractureMatrix{pEDFM_alpha_Tx_index(t)},' ');
                            i     = str2double(pEDFM_alpha_Tx_Split{2})+1;
                            j     = str2double(pEDFM_alpha_Tx_Split{3})+1;
                            k     = str2double(pEDFM_alpha_Tx_Split{4})+1;
                            alpha = str2double(pEDFM_alpha_Tx_Split{5});
                            pEDFM_alpha_Tx(i,j,k) = alpha;
                        end
                        % pEDFM_alpha_Ty
                        temp = strfind(FractureMatrix, 'ROCK_ALPHA_TY');
                        pEDFM_alpha_Ty_index = find(~cellfun('isempty', temp));
                        for t = 1 : length(pEDFM_alpha_Ty_index)
                            pEDFM_alpha_Ty_Split = strsplit(FractureMatrix{pEDFM_alpha_Ty_index(t)},' ');
                            i     = str2double(pEDFM_alpha_Ty_Split{2})+1;
                            j     = str2double(pEDFM_alpha_Ty_Split{3})+1;
                            k     = str2double(pEDFM_alpha_Ty_Split{4})+1;
                            alpha = str2double(pEDFM_alpha_Ty_Split{5});
                            pEDFM_alpha_Ty(i,j,k) = alpha;
                        end
                        % pEDFM_alpha_Tz
                        temp = strfind(FractureMatrix, 'ROCK_ALPHA_TZ');
                        pEDFM_alpha_Tz_index = find(~cellfun('isempty', temp));
                        for t = 1 : length(pEDFM_alpha_Tz_index)
                            pEDFM_alpha_Tz_Split = strsplit(FractureMatrix{pEDFM_alpha_Tz_index(t)},' ');
                            i     = str2double(pEDFM_alpha_Tz_Split{2})+1;
                            j     = str2double(pEDFM_alpha_Tz_Split{3})+1;
                            k     = str2double(pEDFM_alpha_Tz_Split{4})+1;
                            alpha = str2double(pEDFM_alpha_Tz_Split{5});
                            pEDFM_alpha_Tz(i,j,k) = alpha;
                        end
                        ReservoirGrid.AddpEDFMCorrections( pEDFM_alpha_Tx , pEDFM_alpha_Ty , pEDFM_alpha_Tz );
                        
                    case('CornerPointGrid')
                        pEDFM_alpha_T = zeros( length(ReservoirGrid.Trans) , 1 );
                        temp = strfind(FractureMatrix, 'ROCK_pEDFM_ALPHA');
                        pEDFM_alpha_T_index = find(~cellfun('isempty', temp));
                        for t = 1 : length(pEDFM_alpha_T_index)
                            pEDFM_alpha_T_Split = strsplit(FractureMatrix{pEDFM_alpha_T_index(t)},' ');
                            i     = str2double(pEDFM_alpha_T_Split{2})+1;
                            alpha = str2double(pEDFM_alpha_T_Split{3});
                            pEDFM_alpha_T(i) = alpha;
                        end
                        ReservoirGrid.AddpEDFMCorrections( pEDFM_alpha_T );
                end
                % 2. Fractures
                pEDFM_alpha_Tx = cell(FracturesGrid.Nfrac,1);
                pEDFM_alpha_Ty = cell(FracturesGrid.Nfrac,1);
                pEDFM_alpha_Tz = 0;
                for f = 1 : FracturesGrid.Nfrac
                    pEDFM_alpha_Tx{f} = zeros(FracturesGrid.Grids(f).Nx+1,FracturesGrid.Grids(f).Ny);
                    pEDFM_alpha_Ty{f} = zeros(FracturesGrid.Grids(f).Nx,FracturesGrid.Grids(f).Ny+1);
                end
                % Tx_alpha
                temp = strfind(FractureMatrix, 'FRAC_ALPHA_TX');
                pEDFM_alpha_Tx_index = find(~cellfun('isempty', temp));
                for t = 1 : length(pEDFM_alpha_Tx_index)
                    pEDFM_alpha_Tx_Split = strsplit(FractureMatrix{pEDFM_alpha_Tx_index(t)},' ');
                    f   = str2double(pEDFM_alpha_Tx_Split{2})+1;
                    i_f = str2double(pEDFM_alpha_Tx_Split{3})+1;
                    j_f = str2double(pEDFM_alpha_Tx_Split{4})+1;
                    alpha = str2double(pEDFM_alpha_Tx_Split{5});
                    pEDFM_alpha_Tx{f}(i_f,j_f) = alpha;
                end
                % Ty_alpha
                temp = strfind(FractureMatrix, 'FRAC_ALPHA_TY');
                pEDFM_alpha_Ty_index = find(~cellfun('isempty', temp));
                for t = 1 : length(pEDFM_alpha_Ty_index)
                    pEDFM_alpha_Ty_Split = strsplit(FractureMatrix{pEDFM_alpha_Ty_index(t)},' ');
                    f   = str2double(pEDFM_alpha_Ty_Split{2})+1;
                    i_f = str2double(pEDFM_alpha_Ty_Split{3})+1;
                    j_f = str2double(pEDFM_alpha_Ty_Split{4})+1;
                    alpha = str2double(pEDFM_alpha_Ty_Split{5});
                    pEDFM_alpha_Ty{f}(i_f,j_f) = alpha;
                end
                for f = 1 : FracturesGrid.Nfrac
                    FracturesGrid.Grids(f).AddpEDFMCorrections(pEDFM_alpha_Tx{f},pEDFM_alpha_Ty{f},pEDFM_alpha_Tz)
                end
            end
            
            fprintf(' ---> Complete!\n');
        end
    end
end