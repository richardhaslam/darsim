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
            
%             S1=  reshape(VarValues(:, 2),99,99,1);
%             S1(1:50,:) = ones(50,99);
%             S1 = reshape(S1,N,1);
%             VarValues(:, 2) = S1;
% 
%             S2 =  reshape(VarValues(:, 3),99,99,1);
%             S2(51:end,:) = zeros(49,99);
%             S2 = reshape(S2,N,1);
%             VarValues(:, 3) = S2;            
            
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
                case("Geothermal_2T")
                    VarNames = {'P_1', 'Tf', 'Tr' ,'S_1'};
                    VarValues(:, 4) = 1;
                    VarValues(:, 2:3) = obj.SimulationInput.ReservoirProperties.Temperature;
                    simulation.Initializer = initializer_singlephase(VarNames, VarValues);
                case("Geothermal_1T")
                    VarNames = {'P_1', 'T', 'S_1'};
                    VarValues(:, 3) = 1;
                    VarValues(:, 2) = obj.SimulationInput.ReservoirProperties.Temperature;
                    simulation.Initializer = initializer_singlephase(VarNames, VarValues);
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
                temp = strfind(FractureMatrix, 'TYPE');
                fracGen_Type = find(~cellfun('isempty', temp));
                fracGen_Type = strsplit(FractureMatrix{fracGen_Type},' ');
                
                fprintf('This simulation uses input from "%s" fracture generator.\n', fracGen_Type{2});
                NrOfFrac = obj.SimulationInput.FracturesProperties.NrOfFrac;
                fprintf('Extracting data from %02d fractures ...\n', NrOfFrac);
                FracturesGrid = fractures_grid(NrOfFrac);
                
                temp = strfind(FractureMatrix, 'DIMENSION');
                frac_input_res_dimen = find(~cellfun('isempty', temp));
                frac_input_res_dimen = strsplit(FractureMatrix{frac_input_res_dimen},' ');
                if ( obj.SimulationInput.ReservoirProperties.size(1) ~= str2double(frac_input_res_dimen{2}) ) || ...
                        ( obj.SimulationInput.ReservoirProperties.size(2) ~= str2double(frac_input_res_dimen{4}) ) || ...
                        ( obj.SimulationInput.ReservoirProperties.size(3) ~= str2double(frac_input_res_dimen{6}) )
                    error('The dimension of reservoir in the "fracture" input file does not match the simulation input file!');
                end
                
                temp = strfind(FractureMatrix, 'RESERVOIR_GRID');
                frac_input_res_grid = find(~cellfun('isempty', temp));
                frac_input_res_grid = strsplit(FractureMatrix{frac_input_res_grid},' ');
                if ( ReservoirGrid.Nx ~= str2double(frac_input_res_grid{2}) ) || ...
                        ( ReservoirGrid.Ny ~= str2double(frac_input_res_grid{4}) ) || ...
                        ( ReservoirGrid.Nz ~= str2double(frac_input_res_grid{6}) )
                    error('The number of grid cells for reservoir in the "fracture" input file does not match the simulation input file!');
                end
                
                
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
                
                %n_phases = str2double(inputMatrix(obj.Comp_Type + 3)); % Number of phases (useful to define size of some objects)
                n_phases = obj.SimulationInput.FluidProperties.NofPhases;
                fprintf('---> Fracture ');
                
                for f = 1 : NrOfFrac
                    if (f>1),  fprintf(repmat('\b', 1, 9+27));  end
                    fprintf('%04d/%04d',f,NrOfFrac);
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
                        if strcmp(fracGen_Type{2},'pEDFM')
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
                        if strcmp(fracGen_Type{2},'pEDFM')
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
                if strcmp(fracGen_Type{2},'pEDFM')
                    
                    % reservoir
                    Tx_alpha = zeros(ReservoirGrid.Nx+1,ReservoirGrid.Ny,ReservoirGrid.Nz);
                    Ty_alpha = zeros(ReservoirGrid.Nx,ReservoirGrid.Ny+1,ReservoirGrid.Nz);
                    Tz_alpha = zeros(ReservoirGrid.Nx,ReservoirGrid.Ny,ReservoirGrid.Nz+1);
                    % Tx_alpha
                    temp = strfind(FractureMatrix, 'ROCK_ALPHA_TX');
                    ALPHA_TX_index = find(~cellfun('isempty', temp));
                    for t = 1 : length(ALPHA_TX_index)
                        ALPHA_TX_Split = strsplit(FractureMatrix{ALPHA_TX_index(t)},' ');
                        i = str2double(ALPHA_TX_Split{2})+1;
                        j = str2double(ALPHA_TX_Split{3})+1;
                        k = str2double(ALPHA_TX_Split{4})+1;
                        alpha = str2double(ALPHA_TX_Split{5});
                        Tx_alpha(i,j,k) = alpha;
                    end
                    % Ty_alpha
                    temp = strfind(FractureMatrix, 'ROCK_ALPHA_TY');
                    ALPHA_TY_index = find(~cellfun('isempty', temp));
                    for t = 1 : length(ALPHA_TY_index)
                        ALPHA_TY_Split = strsplit(FractureMatrix{ALPHA_TY_index(t)},' ');
                        i = str2double(ALPHA_TY_Split{2})+1;
                        j = str2double(ALPHA_TY_Split{3})+1;
                        k = str2double(ALPHA_TY_Split{4})+1;
                        alpha = str2double(ALPHA_TY_Split{5});
                        Ty_alpha(i,j,k) = alpha;
                    end
                    % Tz_alpha
                    temp = strfind(FractureMatrix, 'ROCK_ALPHA_TZ');
                    ALPHA_TZ_index = find(~cellfun('isempty', temp));
                    for t = 1 : length(ALPHA_TZ_index)
                        ALPHA_TZ_Split = strsplit(FractureMatrix{ALPHA_TZ_index(t)},' ');
                        i = str2double(ALPHA_TZ_Split{2})+1;
                        j = str2double(ALPHA_TZ_Split{3})+1;
                        k = str2double(ALPHA_TZ_Split{4})+1;
                        alpha = str2double(ALPHA_TZ_Split{5});
                        Tz_alpha(i,j,k) = alpha;
                    end
                    ReservoirGrid.AddpEDFMCorrections(Tx_alpha,Ty_alpha,Tz_alpha)
                    
                    % fractures
                    Tx_alpha = cell(FracturesGrid.Nfrac,1);
                    Ty_alpha = cell(FracturesGrid.Nfrac,1);
                    Tz_alpha = 0;
                    for f = 1 : FracturesGrid.Nfrac
                        Tx_alpha{f} = zeros(FracturesGrid.Grids(f).Nx+1,FracturesGrid.Grids(f).Ny);
                        Ty_alpha{f} = zeros(FracturesGrid.Grids(f).Nx,FracturesGrid.Grids(f).Ny+1);
                    end
                    % Tx_alpha
                    temp = strfind(FractureMatrix, 'FRAC_ALPHA_TX');
                    ALPHA_TX_index = find(~cellfun('isempty', temp));
                    for t = 1 : length(ALPHA_TX_index)
                        ALPHA_TX_Split = strsplit(FractureMatrix{ALPHA_TX_index(t)},' ');
                        f   = str2double(ALPHA_TX_Split{2})+1;
                        i_f = str2double(ALPHA_TX_Split{3})+1;
                        j_f = str2double(ALPHA_TX_Split{4})+1;
                        alpha = str2double(ALPHA_TX_Split{5});
                        Tx_alpha{f}(i_f,j_f) = alpha;
                    end
                    % Ty_alpha
                    temp = strfind(FractureMatrix, 'FRAC_ALPHA_TY');
                    ALPHA_TY_index = find(~cellfun('isempty', temp));
                    for t = 1 : length(ALPHA_TY_index)
                        ALPHA_TY_Split = strsplit(FractureMatrix{ALPHA_TY_index(t)},' ');
                        f   = str2double(ALPHA_TY_Split{2})+1;
                        i_f = str2double(ALPHA_TY_Split{3})+1;
                        j_f = str2double(ALPHA_TY_Split{4})+1;
                        alpha = str2double(ALPHA_TY_Split{5});
                        Ty_alpha{f}(i_f,j_f) = alpha;
                    end
                    for f = 1 : FracturesGrid.Nfrac
                        FracturesGrid.Grids(f).AddpEDFMCorrections(Tx_alpha{f},Ty_alpha{f},Tz_alpha)
                    end
                end
                
                fprintf(' ---> Complete!\n');
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
                                    case {'Geothermal_1T','Geothermal_2T'}
                                        prolongationbuilder.BFUpdater = bf_updater_ms_geothermal();
                                    otherwise
                                        prolongationbuilder.BFUpdater = bf_updater_ms();
                                end
                            else % Fractured
                                switch obj.SimulatorSettings.Formulation
                                    case {'Geothermal_1T','Geothermal_2T'}
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
                    if obj.SimulatorSettings.Formulation == "Geothermal_1T" || obj.SimulatorSettings.Formulation == "Geothermal_2T"
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
                    
                    % a.3 Rock Temperature operator builder
                    if obj.SimulatorSettings.Formulation == "Geothermal_2T"
                        switch (ADMSettings.TrInterpolator)
                            case ('Constant')
                                prolongationbuilder = prolongation_builder_constant(ADMSettings.maxLevel(1));
                            otherwise
                                prolongationbuilder = prolongation_builder_MSRockTemperature(ADMSettings.maxLevel(1), ADMSettings.Coarsening(:,:,1));
                                prolongationbuilder.BFUpdater = bf_updater_ms_geothermal();
                                if strcmp(ADMSettings.TrInterpolator, 'Homogeneous')
                                    prolongationbuilder.BFUpdater.MaxContrast = 1;
                                else
                                    prolongationbuilder.BFUpdater.MaxContrast = ADMSettings.BF_MaxContrast;
                                end
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
                        if obj.SimulatorSettings.Formulation == "Geothermal_2T"
                            prolongationbuilder.BFUpdater = bf_updater_ms_P_geothermal();
                        else
                            prolongationbuilder.BFUpdater = bf_updater_ms();
                        end
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
            cr = obj.SimulationInput.ReservoirProperties.Compressibility;
            nx = obj.SimulationInput.ReservoirProperties.Grid.N(1);
            ny = obj.SimulationInput.ReservoirProperties.Grid.N(2);
            nz = obj.SimulationInput.ReservoirProperties.Grid.N(3);
            Cpr = obj.SimulationInput.ReservoirProperties.SpecificHeat;
            RockDensity = obj.SimulationInput.ReservoirProperties.RockDensity;
            BulkModulus = obj.SimulationInput.ReservoirProperties.BulkModulus;
            ShearModulus = obj.SimulationInput.ReservoirProperties.ShearModulus;
            K = ones(nx*ny*nz, 3);
            switch obj.SimulationInput.ReservoirProperties.PermUnit
                case('m2')
                    PermMultiplier = 1;
                case('D')
                    PermMultiplier = 1e-12;
                case('mD')
                    PermMultiplier = 1e-15;
            end
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
            
            % Adding permeability info to the reservoir
            Reservoir.AddPermeabilityPorosity(K, phi);
            if contains(obj.SimulatorSettings.Formulation,'Geothermal')
                Reservoir.AddConductivity(obj.SimulationInput.ReservoirProperties.RockConductivity,obj.SimulationInput.FluidProperties.FluidConductivity);
            end
            Reservoir.Cr = cr;
            Reservoir.Cpr = Cpr;
            Reservoir.rhoRock = RockDensity;
            Reservoir.bulkMod = BulkModulus;
            Reservoir.shearMod = ShearModulus;
            Reservoir.P0 = obj.SimulationInput.Init(1); % Initial Pressure of the Reservoir
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
                switch (obj.SimulationInput.FluidProperties.FluidModel)
                    case {'Geothermal_1T','Geothermal_2T'}
                        temperature = obj.SimulationInput.WellsInfo.Inj(i).Temperature;
                    otherwise
                        temperature = Tres;
                end
                switch (obj.SimulationInput.WellsInfo.Inj(i).Constraint.name)
                    case('pressure')
                        pressure = obj.SimulationInput.WellsInfo.Inj(i).Constraint.value;
                        % temperature = obj.SimulationInput.WellsInfo.Inj(i).Temperature;
                        Injector = injector_pressure(PI, coord, pressure, temperature, n_phases);
                    case('rate')
                        rate = obj.SimulationInput.WellsInfo.Inj(i).Constraint.value;
                        p_init = obj.SimulationInput.Init(1);
                        rate = rate * Reservoir.TotalPV / (3600 * 24); % convert pv/day to m^3/s
                        Injector = injector_rate(PI, coord, rate, p_init, temperature, n_phases);
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
                    FracturesNetwork.Fractures(f).Temp = Tres;                                          % Temperature of each fracture
                    FracturesNetwork.Fractures(f).Cpr = Cpr;
                    FracturesNetwork.Fractures(f).Rho = RockDensity;
                    
                    Porosity = str2double( frac_info_split{6} );                                        % Porosity  of each fracture
                    Permeability = str2double( frac_info_split{7} );                                    % Permeability of each fracture
                    Kx = ones(FracturesGrid.N(f), 1)*Permeability;
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
                case{'Geothermal_1T','Geothermal_2T'}
                    % build the geothermal fluid model
                    switch (obj.SimulationInput.FluidProperties.FluidModel)
                        case{'Geothermal_1T'}
                            FluidModel = Geothermal_1T_fluid_model();
                        case{'Geothermal_2T'}
                            FluidModel = Geothermal_2T_fluid_model();
                            FluidModel.AveragedTemperature = obj.SimulationInput.FluidProperties.AveragedTemperature;
                    end
                    Phase = therm_comp_phase();
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
                case('Geothermal_1T')
                    Formulation = Geothermal_1T_formulation();
                    obj.NofEq = obj.SimulationInput.FluidProperties.NofPhases + 1;
                case('Geothermal_2T')
                    Formulation = Geothermal_2T_formulation();
                    obj.NofEq = obj.SimulationInput.FluidProperties.NofPhases + 2;
                    Formulation.AveragedTemperature = obj.SimulationInput.FluidProperties.AveragedTemperature;
            end
            Formulation.NofPhases = obj.SimulationInput.FluidProperties.NofPhases;
        end
        function TimeDriver = BuildTimeDriver(obj)
            if obj.SimulatorSettings.LTSPlot == 1
                TimeDriver = LTSPrint_TimeLoop_Driver(obj.SimulatorSettings.reports, obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps);
            % elseif obj.SimulatorSettings.ADT_SEQ == 1 
              %   TimeDriver = LTS_TimeLoop_Driver(obj.SimulatorSettings.reports, obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps);
            else
                 TimeDriver = TimeLoop_Driver(obj.SimulatorSettings.reports, obj.SimulationInput.TotalTime, obj.SimulatorSettings.MaxNumTimeSteps);
            end
            % Construct Coupling
            switch(obj.SimulatorSettings.CouplingType)
                case('FIM')
                    % FIM coupling
                    %%%%FIM settings
                    switch obj.SimulatorSettings.Formulation
                        case {'Geothermal_1T','Geothermal_2T'}
                            NLSolver = NL_Solver_geothermal();
                        otherwise
                            NLSolver = NL_Solver();
                    end
                    NLSolver.SystemBuilder = fim_system_builder();
                    switch obj.SimulatorSettings.DiscretizationModel
                        case ('ADM')
                            % Build a different convergence cheker and a proper LS for ADM
                            switch obj.SimulatorSettings.Formulation
                                case ('Geothermal_1T')
                                    ConvergenceChecker = convergence_checker_ADM_geothermal_1T();
                                case ('Geothermal_2T')
                                    ConvergenceChecker = convergence_checker_ADM_geothermal_2T();
                                    ConvergenceChecker.AveragedTemperature = obj.SimulationInput.FluidProperties.AveragedTemperature;
                                otherwise
                                    ConvergenceChecker = convergence_checker_ADM();
                            end
                            
                            ConvergenceChecker.OperatorsAssembler = operators_assembler_fim(obj.NofEq);
                            NLSolver.LinearSolver = linear_solver_ADM(obj.SimulatorSettings.LinearSolver, 1e-6, 500);
                            NLSolver.LinearSolver.OperatorsAssembler = operators_assembler_fim(obj.NofEq);
                            if obj.SimulatorSettings.Formulation == "Geothermal_2T" && ConvergenceChecker.AveragedTemperature == "On"
                                % In this case, two Tf and Tr eqations with be summed up to represent one average temperature
                                % and for OperatorsAssembler, we have only one equation less (average T instead of Tf and Tr).
                                NLSolver.LinearSolver.OperatorsAssembler = operators_assembler_fim(obj.NofEq-1);
                            end
                            if obj.SimulatorSettings.ADMSettings.DLGR
                                % it will change perm during ADM
                                % simulaiton to use upscaled ones
                                NLSolver.LinearSolver.DLGR = 1;
                            end
                        otherwise
                            switch (obj.SimulatorSettings.Formulation)
                                case('Molar')
                                    ConvergenceChecker = convergence_checker_FS_molar();
                                case ('Geothermal_1T')
                                    ConvergenceChecker = convergence_checker_FS_geothermal_1T();
                                case ('Geothermal_2T')
                                    ConvergenceChecker = convergence_checker_FS_geothermal_2T();
                                    ConvergenceChecker.AveragedTemperature = obj.SimulationInput.FluidProperties.AveragedTemperature;
                                otherwise
                                    ConvergenceChecker = convergence_checker_FS();
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
                        case('Geothermal_1T')
                            ConvergenceChecker.NormCalculator = norm_calculator_geothermal_1T();
                        case('Geothermal_2T')
                            ConvergenceChecker.NormCalculator = norm_calculator_geothermal_2T();
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