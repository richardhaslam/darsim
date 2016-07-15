% Builder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
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
        transport
        plotting
        adm
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
            strategy = find(~cellfun('isempty', temp));
            if strategy ~= 0
                obj.CouplingType = 'FIM';
            else
                temp = strfind(SettingsMatrix{1}, 'SEQUENTIAL'); 
                obj.coupling = find(~cellfun('isempty', temp));
                temp = strfind(SettingsMatrix{1}, 'IMPSAT');
                obj.transport = find(~cellfun('isempty', temp));
            end
            
            temp = strfind(SettingsMatrix{1}, 'ADM');
            obj.adm = find(~cellfun('isempty', temp));
            
            %%%%%%%%%%%%%OPTIONS%%%%%%%%%%%%%%%%
            temp = strfind(inputMatrix{1}, 'OUTPUT');
            xv = find(~cellfun('isempty', temp));
            obj.plotting = char(inputMatrix{1}(xv+1)); %Matlab or VTK
            
        end
        function simulation = BuildSimulation(obj, inputMatrix, SettingsMatrix)
            simulation = Reservoir_Simulation();
            simulation.DiscretizationModel = obj.BuildDiscretization(inputMatrix, SettingsMatrix);
            simulation.ProductionSystem = obj.BuildProductionSystem(inputMatrix, simulation.DiscretizationModel);            
            simulation.FluidModel = obj.BuildFluidModel(inputMatrix);
            simulation.Formulation = obj.BuildFormulation(inputMatrix);
            simulation.TimeDriver = obj.BuildTimeDriver(inputMatrix, SettingsMatrix);
            simulation.Summary = obj.BuildSummary();
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
                Kx = ones(Grid.N,1)*10^(-12);
                Ky = ones(Grid.N,1)*10^(-12);
            end
            K = [Kx, Ky];
            Reservoir.AddPermeabilityPorosity(K, phi);
            ProductionSystem.AddReservoir(Reservoir);
            
            % Wells
            Wells = wells();
            %Injectors
            for i=1:length(obj.inj)
                i_init = str2double(inputMatrix(obj.inj(i) + 1));
                i_final = str2double(inputMatrix(obj.inj(i) + 2));
                j_init = str2double(inputMatrix(obj.inj(i) + 3));
                j_final = str2double(inputMatrix(obj.inj(i) + 4));
                coord = [i_init, i_final; j_init, j_final];
                PI = 2000;
                Injector = injector_pressure(PI, coord, str2double(inputMatrix(obj.inj(i) + 6)));
                Wells.AddInjector(Injector);
            end
            %Producers
            for i=1:length(obj.prod)
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
        function FluidModel = BuildFluidModel(obj, inputMatrix)     
           switch (char(inputMatrix(obj.Comp_Type+1)))
               case ('Immiscible')
                   FluidModel = Immiscible_fluid_model(n_phases);
               case ('BlackOil')
                   FluidModel = BO_fluid_model(n_phases, n_comp);
               case ('Compositional')
                   FluidModel = Comp_fluid_model(n_phases, n_comp);
           end
        end
        function Formulation = BuildFormulation(obj, inputMatrix)
            
        end
        function TimeDriver = BuildTimeDriver(obj, inputMatrix, SettingsMatrix)
        end
        function Summary = BuildSummary(obj)
        end
        function Writer = BuildWriter(obj, simulation)
            
        end
    end
end