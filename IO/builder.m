% Builder Builds all objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 8 September 2016
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
            temp = strfind(inputMatrix{1}, 'COMPOSITION TYPE');
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
            %%%%%%%%%%%%%WELLS%%%%%%%%%%%%%%%%
            temp = regexp(inputMatrix{1}, 'INJ\d', 'match');
            obj.inj = find(~cellfun('isempty', temp));
            temp = regexp(inputMatrix{1}, 'PROD\d', 'match');
            obj.prod = find(~cellfun('isempty', temp));
            
            %%%%%%%%%%%%%%%SIMULATOR'S SETTINGS%%%%%%%%%%%
            temp = strfind(SettingsMatrix{1}, 'TIMESTEPS');
            xv = find(~cellfun('isempty', temp));
            obj.MaxNumTimeSteps = str2double(SettingsMatrix{1}(xv+1));
            temp = strfind(SettingsMatrix{1}, 'REPORTS');
            xv = find(~cellfun('isempty', temp));
            obj.reports = str2double(SettingsMatrix{1}(xv+1));
            temp = strfind(SettingsMatrix{1}, 'FIM'); 
            obj.coupling = find(~cellfun('isempty', temp));
            if obj.coupling ~= 0
                obj.CouplingType = 'FIM';
            else
                obj.CouplingType ='Sequential';
                temp = strfind(SettingsMatrix{1}, 'SEQUENTIAL'); 
                obj.coupling = find(~cellfun('isempty', temp));
                temp = strfind(SettingsMatrix{1}, 'IMPSAT');
                obj.transport = find(~cellfun('isempty', temp));
            end
            
            temp = strfind(SettingsMatrix{1}, 'ADM');
            obj.adm = find(~cellfun('isempty', temp));
            temp = strfind(SettingsMatrix{1}, 'FLASH LOOPS');
            obj.flash = find(~cellfun('isempty', temp));
            
            %%%%%%%%%%%%%OPTIONS%%%%%%%%%%%%%%%%
            temp = strfind(SettingsMatrix{1}, 'OUTPUT');
            xv = find(~cellfun('isempty', temp));
            obj.plotting = char(SettingsMatrix{1}(xv+1)); %Matlab or VTK
            
        end
        function simulation = BuildSimulation(obj, inputMatrix, SettingsMatrix)
            simulation = Reservoir_Simulation();
            simulation.DiscretizationModel = obj.BuildDiscretization(inputMatrix, SettingsMatrix);
            simulation.ProductionSystem = obj.BuildProductionSystem(inputMatrix, simulation.DiscretizationModel);            
            simulation.FluidModel = obj.BuildFluidModel(inputMatrix, simulation.ProductionSystem);
            simulation.Formulation = obj.BuildFormulation(inputMatrix, simulation.DiscretizationModel, simulation.FluidModel);
            simulation.TimeDriver = obj.BuildTimeDriver(simulation.FluidModel, SettingsMatrix);
            simulation.Summary = obj.BuildSummary(simulation);
        end
        function Discretization = BuildDiscretization(obj, inputMatrix, SettingsMatrix)
            %Gridding
            nx = str2double(inputMatrix(obj.grid + 1));
            ny = str2double(inputMatrix(obj.grid + 2));
            nz = str2double(inputMatrix(obj.grid + 3));
            if (str2double(SettingsMatrix(obj.adm + 1)) == 0 )
                Discretization = FS_Discretization_model(nx, ny, nz);
                obj.ADM = 'inactive';
            else
                obj.ADM = 'active';
                temp = strfind(SettingsMatrix, 'LEVELS');
                x = find(~cellfun('isempty', temp));
                maxLevel = str2double(SettingsMatrix(x+1));
                temp = strfind(SettingsMatrix, 'TOLERANCE');
                x = find(~cellfun('isempty', temp));
                tol = str2double(SettingsMatrix(x+1));
                temp = strfind(SettingsMatrix, 'COARSENING_RATIOS');
                x = find(~cellfun('isempty', temp));
                cx = str2double(SettingsMatrix(x+1));
                cy = str2double(SettingsMatrix(x+2));
                Coarsening = [cx, cy; cx^2, cy^2; cx^3, cy^3]; %Coarsening Factors: Cx1, Cy1; Cx2, Cy2; ...; Cxn, Cyn;
                temp = strfind(SettingsMatrix, 'PRESSURE_INTERPOLATOR');
                x = find(~cellfun('isempty', temp));
                switch (char(SettingsMatrix(x+1))) 
                    case ('Constant')
                        operatorshandler = operators_handler_constant(maxLevel, cx*cy);
                    case ('Homogeneous')
                        operatorshandler = operators_handler_MS(maxLevel, cx*cy);
                        operatorshandler.BFUpdater = bf_updater_bilin();
                    case ('MS')
                        operatorshandler = operators_handler_MS(maxLevel, cx*cy);
                        operatorshandler.BFUpdater = bf_updater_ms();
                end
                gridselector = adm_grid_selector(tol);
                Discretization = ADM_Discretization_model(nx, ny, nz, maxLevel, Coarsening);
                Discretization.AddADMGridSelector(gridselector);
                Discretization.AddOperatorsHandler(operatorshandler);
            end
            
        end
        function ProductionSystem = BuildProductionSystem (obj, inputMatrix, DiscretizationModel)
            ProductionSystem = Production_System();
            %Reservoir
            Lx = str2double(inputMatrix(obj.size +1));  %Dimension in x−direction [m]
            Ly = str2double(inputMatrix(obj.size +2));  %Dimension in y−direction [m]
            h = str2double(inputMatrix(obj.size + 3));  %Reservoir thickness [m]
            Tres = str2double(inputMatrix(obj.temperature + 1));   %Res temperature [K]
            Reservoir = reservoir(Lx, Ly, h, Tres);
            phi = str2double(inputMatrix(obj.por + 1));
            if strcmp(inputMatrix(obj.perm - 1), 'INCLUDE')
                %File name
                file  = strcat('../Permeability/', char(inputMatrix(obj.perm +1)));
                %load the file in a vector
                field = load(file);
                % reshape it to specified size
                field = reshape(field(3:end),[field(1) field(2)]);
                % make it the size of the grid
                Kx = reshape(field(1:DiscretizationModel.ReservoirGrid.Nx,1:DiscretizationModel.ReservoirGrid.Ny), DiscretizationModel.ReservoirGrid.N, 1);
                Ky = Kx;
            else
                value = str2double(inputMatrix(obj.perm +1));
                Kx = ones(DiscretizationModel.ReservoirGrid.N, 1)*value;
                Ky = ones(DiscretizationModel.ReservoirGrid.N, 1)*value;
            end
            K = [Kx, Ky];
            Reservoir.AddPermeabilityPorosity(K, phi);
            ProductionSystem.AddReservoir(Reservoir);
            
            % Wells
            Wells = wells();
            Wells.NofInj = length(obj.inj);
            Wells.NofProd = length(obj.prod);
            %Injectors
            for i=1:Wells.NofInj
                i_init = str2double(inputMatrix(obj.inj(i) + 1));
                i_final = str2double(inputMatrix(obj.inj(i) + 2));
                j_init = str2double(inputMatrix(obj.inj(i) + 3));
                j_final = str2double(inputMatrix(obj.inj(i) + 4));
                coord = [i_init, i_final; j_init, j_final];
                PI = 2000;
                pressure = str2double(inputMatrix(obj.inj(i) + 6));
                Injector = injector_pressure(PI, coord, pressure, Tres);
                Wells.AddInjector(Injector);
            end
            %Producers
            for i=1:Wells.NofProd
                i_init = str2double(inputMatrix(obj.prod(i) + 1));
                i_final = str2double(inputMatrix(obj.prod(i) + 2));
                j_init = str2double(inputMatrix(obj.prod(i) + 3));
                j_final = str2double(inputMatrix(obj.prod(i) + 4));
                coord = [i_init, i_final; j_init, j_final];
                PI = 2000;
                pressure = str2double(inputMatrix(obj.prod(i) + 6));
                Producer = producer_pressure(PI, coord, pressure);
                Wells.AddProducer(Producer);
            end
            ProductionSystem.AddWells(Wells);
            
        end
        function FluidModel = BuildFluidModel(obj, inputMatrix, ProductionSystem)
            n_phases = str2double(inputMatrix(obj.Comp_Type + 3));
            n_comp = str2double(inputMatrix(obj.Comp_Type + 5));
            switch(char(inputMatrix(obj.Comp_Type+1)))
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
                case('BlackOil')
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
                    FluidModel.KvaluesCalculator = BO_Kvalues_calculator();
                case('Compositional')
                    FluidModel = Comp_fluid_model(n_phases, n_comp);
                    FluidModel.KvaluesCalculator = Constant_Kvalues_calculator();
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
                        %Gets all atmospheric bubble points [K]
                        Tb = str2double(inputMatrix(obj.Comp_Prop + 3 + (i-1)*7));
                        %Gets all slopes connecting bubble point and
                        %critical point on 1/T plot [K]
                        b = str2double(inputMatrix(obj.Comp_Prop + 5 + (i-1)*7));
                        MM = str2double(inputMatrix(obj.Comp_Prop + 7*i));
                        comp = component();
                        comp.AddCompProperties(Tb, b, MM);
                        FluidModel.AddComponent(comp, i); 
                    end
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
        function Formulation = BuildFormulation(obj, inputMatrix, Discretization, FluidModel)
            formulationtype = 'Natural';
            %formulationtype ='OverallComposition';
            if (strcmp(char(inputMatrix(obj.Comp_Type+1)), 'Immiscible') == 1)
                formulationtype = 'Immiscible';
            end
            if (strcmp(obj.CouplingType, 'Sequential') == 1)
                formulationtype = 'Sequential';
            end
            switch(formulationtype)
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
                case('OverallComposition')
                    Formulation = Overall_Composition_formulation(FluidModel.NofComp);
                    if strcmp(obj.ADM, 'active')
                        Discretization.OperatorsHandler.FullOperatorsAssembler = operators_assembler_Imm();
                    end
                case('Sequential')
                    Formulation = seq_formulation();
            end
            
            % Gravity model
            Formulation.GravityModel = gravity_model(Discretization.ReservoirGrid.Nx, Discretization.ReservoirGrid.Ny);
            switch (obj.Gravity)
                case('ON')
                    Formulation.GravityModel.g = 9.806;
                case('OFF')
                    Formulation.GravityModel.g = 0;
            end
        end
        function TimeDriver = BuildTimeDriver(obj, FluidModel, SettingsMatrix)
            TimeDriver = TimeLoop_Driver(obj.TotalTime, obj.reports);
            TimeDriver.MaxNumberOfTimeSteps = obj.MaxNumTimeSteps;
            %% Construct Coupling
            switch(obj.CouplingType)
                case('FIM')
                    %%%%FIM settings
                    % Build a different convergence cheker and a proper LS for ADM
                     NLSolver = NL_Solver();
                    if (str2double(SettingsMatrix(obj.adm + 1))==0)
                        NLSolver.SystemBuilder = fim_system_builder();
                        ConvergenceChecker = convergence_checker_FS();
                        NLSolver.LinearSolver = linear_solver();
                    else
                        NLSolver.SystemBuilder = fim_system_builder_ADM();
                        ConvergenceChecker = convergence_checker_ADM();
                        NLSolver.LinearSolver = linear_solver_ADM();
                    end
                    NLSolver.MaxIter = str2double(SettingsMatrix(obj.coupling + 1));
                    ConvergenceChecker.Tol = str2double(SettingsMatrix(obj.coupling + 2));
                    switch (FluidModel.name)
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
                    pressuresolver = incompressible_pressure_solver();
                    pressuresolver.LinearSolver = linear_solver();
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
                        transportsolver.LinearSolver = linear_solver();
                    else
                        transportsolver = explicit_transport_solver();
                        transportsolver.SystemBuilder = explicit_transport_system_builder();
                        Coupling.ConvergenceChecker = convergence_checker_impes();
                    end
                    Coupling.AddTransportSolver(transportsolver);
            end
            Coupling.TimeStepSelector = timestep_selector(str2double(SettingsMatrix(obj.coupling + 3)));
            TimeDriver.AddCouplingStrategy(Coupling);
        end
        function Summary = BuildSummary(obj, simulation)
            %%%%%%%%%%%%%%% BuildObjects for OUTPUT%%%%%%%%%
            switch(obj.CouplingType)
                case('FIM')
                    CouplingStats = FIM_Stats(obj.MaxNumTimeSteps);
                case('Sequential')
                    CouplingStats = Sequential_Stats(obj.MaxNumTimeSteps);
            end
            wellsData = wells_data(obj.MaxNumTimeSteps, simulation.FluidModel.NofPhases, simulation.FluidModel.NofComp, simulation.ProductionSystem.Wells);
            Summary = Run_Summary(obj.MaxNumTimeSteps, CouplingStats, wellsData);
        end
        function Writer = BuildWriter(obj, InputDirectory, simulation)
            % Build Plotter
            switch(obj.plotting)
                case('Matlab')
                    if simulation.DiscretizationModel.ReservoirGrid.Nx == 1 || simulation.DiscretizationModel.ReservoirGrid.Ny == 1
                        plotter = Matlab_Plotter_1D(simulation.ProductionSystem.Wells.Prod.p, simulation.ProductionSystem.Wells.Inj.p);
                    else
                        plotter = Matlab_Plotter_2D(simulation.ProductionSystem.Wells.Prod(1).p, simulation.ProductionSystem.Wells.Inj(1).p);
                    end
                case('VTK')
                    plotter = VTK_Plotter(InputDirectory, obj.ProblemName);
                otherwise
                    plotter = no_Plotter();
            end
            
            switch(obj.ADM)
                case ('inactive')
                    Writer = output_writer_FS(InputDirectory, obj.ProblemName, simulation.ProductionSystem.Wells.NofInj, simulation.ProductionSystem.Wells.NofProd, simulation.Summary.CouplingStats.NTimers, simulation.Summary.CouplingStats.NStats);
                case ('active')
                    Writer = output_writer_adm(InputDirectory, obj.ProblemName, simulation.ProductionSystem.Wells.NofInj, simulation.ProductionSystem.Wells.NofProd, simulation.Summary.CouplingStats.NTimers, simulation.Summary.CouplingStats.NStats);
            end
            Writer.AddPlotter(plotter);
        end
    end
end