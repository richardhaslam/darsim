% Builder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 18 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef builder < handle
    properties
        ProblemName
        TotalTime
        size
        temperature
        grid
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
        CouplingType
        coupling
        transport
        plotting
        adm
        flash
    end
    methods
        function FindKeyWords(obj, inputMatrix, SettingsMatrix)
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
            temp = strfind(inputMatrix{1}, 'DENSITY (kg/m^3)'); 
            obj.density = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'VISCOSITY (Pa sec)');
            obj.viscosity = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'RELPERM');
            obj.relperm = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'COMPRESSIBILITY (1/Pa)');
            obj.compressibility = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'CAPILLARITY');
            obj.capillarity = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'FOAM');
            obj.foam = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'COMPOSITION TYPE');
            obj.Comp_Type = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'COMPONENT PROPERTIES');
            obj.Comp_Prop = find(~cellfun('isempty', temp));
            
            %%%%%%%%%%%%%WELLS%%%%%%%%%%%%%%%%
            temp = regexp(inputMatrix{1}, 'INJ\d', 'match');
            obj.inj = find(~cellfun('isempty', temp));
            temp = regexp(inputMatrix{1}, 'PROD\d', 'match');
            obj.prod = find(~cellfun('isempty', temp));
            
            %%%%%%%%%%%%%%%SIMULATOR'S SETTINGS%%%%%%%%%%%
            temp = strfind(SettingsMatrix{1}, 'TIMESTEPS');
            xv = find(~cellfun('isempty', temp));
            obj.MaxNumTimeSteps = str2double(SettingsMatrix{1}(xv+1));
            temp = strfind(SettingsMatrix{1}, 'FIM'); 
            obj.coupling = find(~cellfun('isempty', temp));
            if obj.coupling ~= 0
                obj.CouplingType = 'FIM';
            else
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
            simulation.FluidModel = obj.BuildFluidModel(inputMatrix, SettingsMatrix);
            simulation.Formulation = obj.BuildFormulation(inputMatrix);
            simulation.TimeDriver = obj.BuildTimeDriver(SettingsMatrix);
            simulation.Summary = obj.BuildSummary(simulation);
        end
        function Discretization = BuildDiscretization(obj, inputMatrix, SettingsMatrix)
            %Gridding
            nx = str2double(inputMatrix(obj.grid + 1));
            ny = str2double(inputMatrix(obj.grid + 2));
            nz = str2double(inputMatrix(obj.grid + 3));
            if (str2double(SettingsMatrix(obj.adm + 1)) == 0 )
                Discretization = FS_Discretization_model(nx, ny, nz);
            else
                temp = strfind(inputMatrix{1}, 'PRESSURE_INTERPOLATOR');
                x = find(~cellfun('isempty', temp));
                ADMSettings.Pressure_Interpolator =  char(inputMatrix{1}(x+1));
                temp = strfind(inputMatrix{1}, 'LEVELS');
                x = find(~cellfun('isempty', temp));
                ADMSettings.maxLevel = str2double(inputMatrix{1}(x+1));
                temp = strfind(inputMatrix{1}, 'TOLERANCE');
                x = find(~cellfun('isempty', temp));
                ADMSettings.tol = str2double(inputMatrix{1}(x+1));
                temp = strfind(inputMatrix{1}, 'COARSENING_RATIOS');
                x = find(~cellfun('isempty', temp));
                cx = str2double(inputMatrix{1}(x+1));
                cy = str2double(inputMatrix{1}(x+2));
                ADMSettings.Coarsening = [cx, cy; cx^2, cy^2; cx^3, cy^3]; %Coarsening Factors: Cx1, Cy1; Cx2, Cy2; ...; Cxn, Cyn;
                Discretization = ADM_Discretization_model(nx, ny, nz, ADMSettings);
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
                Kx = reshape(field(1:DiscretizationModel.ReservoirGrid.Nx,1:DiscretizationModel.ReservoirGrid.Ny)*10^(-12), DiscretizationModel.ReservoirGrid.N, 1);
                Ky = Kx;
            else
                Kx = ones(DiscretizationModel.ReservoirGrid.N,1)*10^(-12);
                Ky = ones(DiscretizationModel.ReservoirGrid.N,1)*10^(-12);
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
                Producer = producer_pressure(PI, coord, str2double(inputMatrix(obj.prod(i) + 6)));
                Wells.AddProducer(Producer);
            end
            ProductionSystem.AddWells(Wells);
            
        end
        function FluidModel = BuildFluidModel(obj, inputMatrix, SettingsMatrix)
            n_phases = str2double(inputMatrix(obj.Comp_Type + 3));
            n_comp = str2double(inputMatrix(obj.Comp_Type + 5));
            switch(char(inputMatrix(obj.Comp_Type+1)))
                case('Immiscible')
                    FluidModel = Immiscible_fluid_model(n_phases);
                case('BlackOil')
                    FluidModel = BO_fluid_model(n_phases, n_comp);
                    FluidModel.Pref = 1e5;
                case('Compositional')
                    FluidModel = Comp_fluid_model(n_phases, n_comp);
                    % Add components
                    for i = 1:FluidModel.NofComp
                        %Gets all atmospheric bubble points [K]
                        Tb = str2double(inputMatrix(obj.Comp_Prop + 3 + (i-1)*5));
                        %Gets all slopes connecting bubble point and
                        %critical point on 1/T plot [K]
                        b = str2double(inputMatrix(obj.Comp_Prop + 5*i));
                        comp = component();
                        comp.AddCompProperties(Tb, b);
                        FluidModel.AddComponent(comp, i); 
                    end
                                        
                    FlashSettings.TolInner = str2double(SettingsMatrix(obj.flash + 2));
                    FlashSettings.MaxIt = str2double(SettingsMatrix(obj.flash + 3));
                    FlashSettings.TolFlash = str2double(SettingsMatrix(obj.flash + 4));
                    FluidModel.AddFlash(FlashSettings);
            end
            % Add phases
            for i = 1:FluidModel.NofPhases
                Phase = phase();
                %Gets all densities [kg/m^3]
                Phase.rho0 = str2double(inputMatrix(obj.density + 2*i));
                %Gets all viscosities [Pa sec]
                Phase.mu = str2double(inputMatrix(obj.viscosity + 2*i));
                %Gets all compressibilities [1/Pa]
                Phase.cf = str2double(inputMatrix(obj.compressibility + 2*i));
                %Gets all residual saturations [-]
                Phase.sr = str2double(inputMatrix(obj.relperm + 1 + 2*i));
                FluidModel.AddPhase(Phase, i);
            end
            
            
            switch(char(inputMatrix(obj.relperm + 1)))
                case('Linear')
                    FluidModel.RelPermModel = relperm_model_linear();
                case('Quadratic')
                    FluidModel.RelPermModel = relperm_model_quadratic();
            end
            switch(char(inputMatrix(obj.capillarity + 1)))
                case('JLeverett')
                    FluidModel.CapillaryModel = capillary_model_leverett();
                case('BrooksCorey')
                    FluidModel.CapillaryModel = capillary_model_corey();
                case('Linear')
                    FluidModel.CapillaryModel = capillary_model_linear();
            end
            
        end
        function Formulation = BuildFormulation(obj, inputMatrix)
            formulationtype = 'Natural';
            if (strcmp(char(inputMatrix(obj.Comp_Type+1)), 'Immiscible') == 1)
                formulationtype = 'Immiscible';
            end
            switch(formulationtype)
                case('Immiscible')
                    Formulation = Immiscible_formulation();
                case('Natural')
                    Formulation = NaturalVar_formulation();
                case('Mass')
                    Formulation = Mass_formulation();
            end
        end
        function TimeDriver = BuildTimeDriver(obj, SettingsMatrix)
            n_reports = 10; %Hard coded for now
            TimeDriver = TimeLoop_Driver(obj.TotalTime, n_reports);
            TimeDriver.MaxNumberOfTimeSteps = obj.MaxNumTimeSteps;
            %% Construct Coupling
            switch(obj.CouplingType)
                case('FIM')
                    %%%%FIM settings
                    NLSolver = NL_Solver();
                    NLSolver.MaxIter = str2double(SettingsMatrix(obj.coupling + 1));
                    % Build a different convergence cheker and a proper LS for ADM
                    if (str2double(SettingsMatrix(obj.adm + 1))==0)
                        ConvergenceChecker = convergence_checker_FS();
                        NLSolver.LinearSolver = linear_solver();
                    else
                        ConvergenceChecker = convergence_checker_ADM();
                        NLSolver.LinearSolver = linear_solver_ADM();
                    end
                    ConvergenceChecker.Tol = str2double(SettingsMatrix(obj.coupling + 2));                 
                    NLSolver.AddConvergenceChecker(ConvergenceChecker); 
                    
                    Coupling = FIM_Strategy('FIM', NLSolver);
                    Coupling.CFL = str2double(SettingsMatrix(obj.coupling + 3));  
                case('Sequential')
                    Coupling = Sequential_Strategy();
            end
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
                        plotter = Matlab_Plotter_1D();
                    else
                        plotter = Matlab_Plotter_2D();
                    end
                case('VTK')
                    plotter = VTK_Plotter(InputDirectory, obj.ProblemName);
                otherwise
                    plotter = no_Plotter();
            end
            
            Writer = output_writer_txt(InputDirectory, obj.ProblemName, simulation.ProductionSystem.Wells.NofInj, simulation.ProductionSystem.Wells.NofProd, simulation.Summary.CouplingStats.NTimers, simulation.Summary.CouplingStats.NStats);
            Writer.AddPlotter(plotter);
        end
    end
end