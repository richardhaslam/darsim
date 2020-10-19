% Geothermal formulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Geothermal_SinglePhase_formulation < formulation
    properties
        MatrixAssembler
        Thph % Transmisibility matrix for mas convcection flux
        Ghph % Gravity         matrix for mas convcection flux
        Tk   % Transmisibility matrix for heat conduction flux
        drhodT
        d2rhodT2
        dmudT
        d2mudT2
        dhdp
        dhdT
        d2hdT2
        dMobdT
        d2MobdT2
        dMobdp
        Cp_std % fluid specific heat in standard condition
        Kf % fluid conductivity
    end
    methods
        function obj = Geothermal_SinglePhase_formulation()
            obj@formulation();
            obj.Tph = cell(1,1);
            obj.Gph = cell(1,1);
            obj.Tk = cell(1,1);
            obj.MatrixAssembler = matrix_assembler_geothermal();
        end
        function x = GetPrimaryUnknowns(obj, ProductionSystem, DiscretizationModel)
             Nt = DiscretizationModel.N; % total grid
             Nm = DiscretizationModel.ReservoirGrid.N; % matrix grid
             Nf = Nt - Nm; % fracture grid
             x = zeros(3 * Nt, 1); % 3 unknowns: P_1, T, Tr
             Index.Start = 1;
             Index.End = Nm;
                     
             for ph=1:obj.NofPhases
                 % Get matrix unknowns
                 x(Index.Start     :Index.End       ) = ProductionSystem.Reservoir.State.Properties(['P_',num2str(ph)]).Value;
                 x(Nt+Index.Start  :Nt + Index.End  ) = ProductionSystem.Reservoir.State.Properties('T').Value;
                 
                 % Get fracture unknowns
                 for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                     Index.Start = Index.End+1;
                     Index.End = Index.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                     x(Index.Start     :Index.End       ) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_',num2str(ph)]).Value;
                     x(Nt+Index.Start  :Nt + Index.End  ) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('T').Value;   
                 end
             end
        end
        function x = GetPrimaryPressure(obj, ProductionSystem, DiscretizationModel)
            Nt = DiscretizationModel.N;
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nf = Nt - Nm;
            x = zeros(Nt, 1);
            if obj.NofPhases > 1
                x(1:Nm) = ProductionSystem.Reservoir.State.Properties('P_2').Value;
            else
                x(1:Nm) = ProductionSystem.Reservoir.State.Properties('P_1').Value;
            end
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            obj.ComputeProperties(ProductionSystem, FluidModel);
            obj.ComputeDerivatives(ProductionSystem, FluidModel);
        end
        function ComputeProperties(obj, ProductionSystem, FluidModel)
            %% 0. Geothermal Properties
            obj.Cp_std = FluidModel.Phases.Cp_std; % define Cp from simualtion builder
            obj.Kf     = FluidModel.Phases.Kf; % define kf from simualtion builder
            %% 1. Reservoir Properties
            FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
            FluidModel.ComputePhaseViscosities(ProductionSystem.Reservoir.State);
            FluidModel.ComputePhaseEnthalpies(ProductionSystem.Reservoir.State);
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State);
            FluidModel.ComputeThermalConductivity(ProductionSystem.Reservoir)
            %% 2. Fractures Properties
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                FluidModel.ComputePhaseViscosities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                FluidModel.ComputePhaseEnthalpies(ProductionSystem.FracturesNetwork.Fractures(f).State);
                obj.Mob = [obj.Mob; FluidModel.ComputePhaseMobilities(ProductionSystem.FracturesNetwork.Fractures(f).State)];
                FluidModel.ComputeThermalConductivity(ProductionSystem.FracturesNetwork.Fractures(f))
            end
        end
        function ComputeDerivatives(obj, ProductionSystem, FluidModel)
            %% 1. Reservoir Derivatives
            obj.drhodp = FluidModel.ComputeDrhoDp(ProductionSystem.Reservoir.State);
            [obj.drhodT,obj.d2rhodT2] = FluidModel.ComputeDrhoDT(ProductionSystem.Reservoir.State);
            [obj.dmudT,obj.d2mudT2] = FluidModel.ComputeDmuDT(ProductionSystem.Reservoir.State);
            obj.dhdp = FluidModel.ComputeDhDp(ProductionSystem.Reservoir.State);
            [obj.dhdT,obj.d2hdT2] = FluidModel.ComputeDhDT(ProductionSystem.Reservoir.State);
            [obj.dMobdT,obj.d2MobdT2] = FluidModel.ComputeDMobdT(ProductionSystem.Reservoir.State);
            obj.dMobdp = zeros(size(obj.drhodp));
             
            %% 2. Fractures Derivatives
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                obj.drhodp = [obj.drhodp; FluidModel.ComputeDrhoDp(ProductionSystem.FracturesNetwork.Fractures(f).State)];
                [drhodT_f,d2rhodT2_f] = FluidModel.ComputeDrhoDT(ProductionSystem.FracturesNetwork.Fractures(f).State);
                obj.drhodT = [obj.drhodT; drhodT_f];
                obj.d2rhodT2 = [obj.d2rhodT2; d2rhodT2_f];
                [dmudT_f, d2mudT2_f] = FluidModel.ComputeDmuDT(ProductionSystem.FracturesNetwork.Fractures(f).State);
                obj.dmudT = [obj.dmudT; dmudT_f];
                obj.d2mudT2 = [obj.d2mudT2; d2mudT2_f];
                obj.dhdp = [obj.dhdp; FluidModel.ComputeDhDp(ProductionSystem.FracturesNetwork.Fractures(f).State)];
                [dhdT_f, d2hdT2_f] = FluidModel.ComputeDhDT(ProductionSystem.FracturesNetwork.Fractures(f).State);
                obj.dhdT = [obj.dhdT; dhdT_f];
                obj.d2hdT2 = [obj.d2hdT2; d2hdT2_f];
                [dMobdT_f,d2MobdT2_f] = FluidModel.ComputeDMobdT(ProductionSystem.FracturesNetwork.Fractures(f).State);
                obj.dMobdT = [obj.dMobdT; dMobdT_f];
                obj.d2MobdT2 = [obj.d2MobdT2; d2MobdT2_f];
            end
        end
        function [Residual_MB   , RHS_MB  ] = BuildMediumMassBalanceResidual(obj, Medium, Grid, dt, State_old, Index, WellMassFlux, FracturesMassFlux, f)
            % Variables at previous time-step n
            P_old   = State_old.Properties('P_1'  ).Value(Index.Start:Index.End);
            rho_old = State_old.Properties('rho_1').Value(Index.Start:Index.End);
            Medium.ComputePorosity(P_old);
            PoreVolume_old = Medium.Por .* Grid.Volume;
            
            % Variables at current time-step iteration stage \nu
            P_new   = Medium.State.Properties('P_1'  ).Value;
            rho_new = Medium.State.Properties('rho_1').Value;
            Medium.ComputePorosity(P_new);
            PoreVolume_new = Medium.Por .* Grid.Volume;
            
            % Mass Accumulation Term
            MassAccumulation = ( PoreVolume_new.*rho_new - PoreVolume_old.*rho_old ) / dt;

            % Mass Balance Residual Calculation (LHS - RHS)
            RHS_MB = WellMassFlux( Index.Start:Index.End , 1 );
            Residual_MB = MassAccumulation ...
                + obj.Tph{1,1+f} * P_new ...
                - obj.Gph{1,1+f} * Grid.Depth ...
                - WellMassFlux( Index.Start:Index.End , 1 ) ...
                - FracturesMassFlux( Index.Start:Index.End , 1 );
        end
        function [Residual_EB   , RHS_EB  ] = BuildMediumEnergyBalanceResidual(obj, Medium, Grid, dt, State_old, Index, WellHeatConvectionFlux, FracturesHeatConvectionFlux, FracturesHeatConductionFlux, f)
            % Variables at previous time-step n
            P_old   = State_old.Properties('P_1'  ).Value(Index.Start:Index.End);
            rho_old = State_old.Properties('rho_1').Value(Index.Start:Index.End);
            T_old   = State_old.Properties('T'    ).Value(Index.Start:Index.End);
            Medium.ComputePorosity(P_old);
            PoreVolume_old = Medium.Por .* Grid.Volume;
            RockVolume_old = (1 - Medium.Por) .* Grid.Volume;
            
            % Variables at current time-step iteration stage \nu
            P_new   = Medium.State.Properties('P_1'  ).Value;
            rho_new = Medium.State.Properties('rho_1').Value;
            T_new   = Medium.State.Properties('T'    ).Value;
            Medium.ComputePorosity(P_new);
            PoreVolume_new = Medium.Por .* Grid.Volume;
            RockVolume_new = (1 - Medium.Por) .* Grid.Volume;
            
            Rho_rock = Medium.Rho;

            % Heat Accumulation Term (U_eff is the effective internal energy)
            U_eff_old = ( rho_old .* obj.Cp_std .* PoreVolume_old + Rho_rock .* Medium.Cpr .* RockVolume_old ) .* T_old;
            U_eff_new = ( rho_new .* obj.Cp_std .* PoreVolume_new + Rho_rock .* Medium.Cpr .* RockVolume_new ) .* T_new;
            HeatAccumulation = ( U_eff_new - U_eff_old ) / dt;
            
            % Energy Balance Residual Calculation (LHS - RHS)
            RHS_EB = WellHeatConvectionFlux( Index.Start:Index.End, 1 );
            Residual_EB  = HeatAccumulation ...
                + obj.Thph{ph, 1+f} * P_new ...
                - obj.Gph{ph, 1+f} * Grid.Depth ...
                + obj.Tk{ph, 1+f} * T_new ...
                - WellHeatConvectionFlux     ( Index.Start:Index.End , 1 )...
                - FracturesHeatConvectionFlux( Index.Start:Index.End , 1 )...
                - FracturesHeatConductionFlux( Index.Start:Index.End , 1 );
        end
        function [Residual_Full , RHS_Full] = BuildFullResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Computing source terms and the flux transfer functions between reservoir and fractures (if any)
            [WellMassFlux, WellHeatConvectionFlux] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            FracturesMassFlux = zeros(DiscretizationModel.N, 1); % Mass flux between each two media (reservoir-frac or frac1-frac2)
            FracturesHeatConvectionFlux = zeros(DiscretizationModel.N, 1); % Heat Convection flux betweem each two media
            FracturesHeatConductionFlux = zeros(DiscretizationModel.N, 1);                           % Heat Conduction flux betweem each two media
            if ProductionSystem.FracturesNetwork.Active
                [FracturesMassFlux, FracturesHeatConvectionFlux, FracturesHeatConductionFlux] = obj.ComputeSourceTerms_frac_mat(ProductionSystem, DiscretizationModel);
            end
            
            % Initializing the Full Residual Vector (N_eq * N_grids, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;        % Number of grid cells in the reservoir
            if ProductionSystem.FracturesNetwork.Active
                Nf = DiscretizationModel.FracturesGrid.N;    % Number of grid cells in the fractures
            else
                Nf = 0;
            end
            Nt = DiscretizationModel.N;                      % Number of grid cells in all media combined
            RHS_Full = zeros( 2*Nt , 1 );
            Residual_Full = zeros( 2*Nt , 1 );
            
            %% Computing Transmissibilities for reservoir
            Index.Start = 1;
            Index.End = Nm;
            for ph=1:obj.NofPhases
                Index.Start = 1;
                Index.End = Nm;
                % Compute mass flow transmissibility and gravity
                [obj.Tph{ph, 1}, obj.Gph{ph, 1}] = obj.MatrixAssembler.TransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
                % Compute heat convection transmissibility and gravity
                [obj.Thph{ph, 1}, obj.Ghph{ph, 1}] = obj.MatrixAssembler.ConvectiveHeatTransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    ProductionSystem.Reservoir.State.Properties(['h_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
            end
            % Compute heat conduction transmissibility
            obj.Tk{1,1} = obj.MatrixAssembler.ConductiveHeatTransmissibilityMatrix( DiscretizationModel.ReservoirGrid , ProductionSystem.Reservoir.State );
            
            %% Computing Transmissibilities for fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Index.Start = Index.End+1;
                Index.End = Index.Start + Nf(f) - 1;
                for ph=1:obj.NofPhases
                    % Compute mass flow transmissibility and gravity
                    [obj.Tph{ph, 1+f}, obj.Gph{ph, 1+f}] = obj.MatrixAssembler.TransmissibilityMatrix( ...
                        DiscretizationModel.FracturesGrid.Grids(f), ...
                        obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph), ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value, ...
                        obj.GravityModel.RhoInt{ph, 1+f} );
                    % Compute heat convection transmissibility and gravity
                    [obj.Thph{ph, 1+f}, obj.Ghph{ph, 1+f}] = obj.MatrixAssembler.ConvectiveHeatTransmissibilityMatrix( ...
                        DiscretizationModel.FracturesGrid.Grids(f), ...
                        obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph), ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value, ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['h_',num2str(ph)]).Value, ...
                        obj.GravityModel.RhoInt{ph, 1+f} );
                end
                % Compute heat conduction transmissibility
                obj.Tk{1,1+f} = obj.MatrixAssembler.ConductiveHeatTransmissibilityMatrix( DiscretizationModel.FracturesGrid.Grids(f) , ProductionSystem.FracturesNetwork.Fractures(f).State );
            end
            
            %% Computing Mass Balance Residual
            % Reservoir
            Index.Start = 1;
            Index.End = Nm;
            [Residual_MB_Reservoir, RHS_MB_Reservoir] = BuildMediumMassBalanceResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, WellMassFlux, FracturesMassFlux, 0);
            Residual_Full( Index.Start : Index.End ) = Residual_MB_Reservoir;
            RHS_Full(      Index.Start : Index.End ) = RHS_MB_Reservoir;
            % Fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Index.Start = Index.End+1;
                Index.End = Index.Start + Nf(f) - 1;
                [Residual_MB_Fractures, RHS_MB_Fractures] = BuildMediumMassBalanceResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, WellMassFlux, FracturesMassFlux, f);
                Residual_Full( Index.Start : Index.End ) = Residual_MB_Fractures;
                RHS_Full(      Index.Start : Index.End ) = RHS_MB_Fractures;
            end
            
            %% Computing Energy Balance Residual
            % Reservoir
            Index.Start = Index.End + 1;
            Index.End = Index.Start + Nm - 1;
            Index_temp.Start = 1;
            Index_temp.End = Nm;
            [Residual_EB_Reservoir, RHS_EB_Reservoir] = BuildMediumEnergyBalanceResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index_temp, WellHeatConvectionFlux, FracturesHeatConvectionFlux, FracturesHeatConductionFlux, 0);
            Residual_Full( Index.Start : Index.End ) = Residual_EB_Reservoir;
            RHS_Full(      Index.Start : Index.End ) = RHS_EB_Reservoir;
            % Fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Index.Start = Index.End + 1;
                Index.End = Index.Start + Nf(f) - 1;
                Index_temp.Start = Index_temp.End + 1;
                Index_temp.End = Index_temp.Start + Nf(f) - 1;
                [Residual_EB_Fractures, RHS_EB_Fractures] = BuildMediumEnergyBalanceResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index_temp, WellHeatConvectionFlux, FracturesHeatConvectionFlux, FracturesHeatConductionFlux, f);
                Residual_Full(Index.Start:Index.End) = Residual_EB_Fractures;
                RHS_Full(     Index.Start:Index.End) = RHS_EB_Fractures;
            end
        end
        function [J_MB_P , J_MB_T] = BuildMediumMassBalanceJacobian(obj, Medium, Wells, Grid, dt, Index, f)
            % Create local variables
            N = Grid.N;
            P = Medium.State.Properties('P_1').Value;
            rho = Medium.State.Properties('rho_1').Value;
            Medium.ComputeDerPorosity(P);
            por = Medium.Por;
            dpor = Medium.DPor;
            
            %% 1.J_PP
            % 1.a Pressure Block
            J_MB_P = obj.Tph{1,1+f};
            switch class(Grid)
                case('corner_point_grid')
                    N_Face = length(obj.UpWind{1,1+f});
                    nc = Grid.N;
                    nf = length(Grid.Trans);
                    C = [ Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index , Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index ];
                    D1 = [ -double(obj.U{1,1+f}>=0)+double(obj.U{1,1+f}<0)  , double(obj.U{1,1+f}>=0)-double(obj.U{1,1+f}<0)]; D1(D1==1)=0;
                    D2 = [ -double(obj.U{1,1+f}>=0)+double(obj.U{1,1+f}<0)  , double(obj.U{1,1+f}>=0)-double(obj.U{1,1+f}<0)];
                    UpwindPermutation1 = sparse([(1:nf)'; (1:nf)'], C, D1, nf, nc)';
                    UpwindPermutation2 = sparse([(1:nf)'; (1:nf)'], C, D2, nf, nc)';
                    
                    % 1.b: compressibility part
                    Mob_drhodp = obj.Mob(Index.Start:Index.End, 1) .* obj.drhodp(Index.Start:Index.End, 1);
                    dMobUpwind = obj.UpWind{1,1+f} * Mob_drhodp;
                    FluxDerivative = dMobUpwind .* abs(obj.U{1,1+f});
                    acc = Grid.Volume/dt .* (por .* obj.drhodp(Index.Start:Index.End) + rho .*dpor);
                    J_MB_P = J_MB_P + spdiags(acc,0,N,N) + (UpwindPermutation1 * spdiags(FluxDerivative,0,N_Face,N_Face) * UpwindPermutation2')';
                    
                    % 2. J_PT
                    dMobdT_drhodT = obj.dMobdT(Index.Start:Index.End, 1) .* rho + obj.Mob(Index.Start:Index.End, 1) .* obj.drhodT(Index.Start:Index.End);
                    dMobUpwind = obj.UpWind{1,1+f} * dMobdT_drhodT;
                    FluxDerivative = dMobUpwind .* abs(obj.U{1,1+f});
                    acc =  Grid.Volume/dt .* por .* obj.drhodT(Index.Start:Index.End) ;
                    J_MB_T = spdiags(acc,0,N,N) + (UpwindPermutation1 * spdiags(FluxDerivative,0,N_Face,N_Face) * UpwindPermutation2')';
                    
                case('cartesian_grid')
                    Nx = Grid.Nx;
                    Ny = Grid.Ny;
                    Nz = Grid.Nz;
                    % 1.b: compressibility part
                    dMupx = obj.UpWind{1,1+f}.x * ( obj.Mob(Index.Start:Index.End, 1) .* obj.drhodp(Index.Start:Index.End) );
                    dMupy = obj.UpWind{1,1+f}.y * ( obj.Mob(Index.Start:Index.End, 1) .* obj.drhodp(Index.Start:Index.End) );
                    dMupz = obj.UpWind{1,1+f}.z * ( obj.Mob(Index.Start:Index.End, 1) .* obj.drhodp(Index.Start:Index.End) );
                    
                    vecX1 = min(reshape(obj.U{1,1+f}.x(1:Nx  ,:     ,:     ), N, 1), 0) .* dMupx;
                    vecX2 = max(reshape(obj.U{1,1+f}.x(2:Nx+1,:     ,:     ), N, 1), 0) .* dMupx;
                    vecY1 = min(reshape(obj.U{1,1+f}.y(:     ,1:Ny  ,:     ), N, 1), 0) .* dMupy;
                    vecY2 = max(reshape(obj.U{1,1+f}.y(:     ,2:Ny+1,:     ), N, 1), 0) .* dMupy;
                    vecZ1 = min(reshape(obj.U{1,1+f}.z(:     ,:     ,1:Nz  ), N, 1), 0) .* dMupz;
                    vecZ2 = max(reshape(obj.U{1,1+f}.z(:     ,:     ,2:Nz+1), N, 1), 0) .* dMupz;
                    acc = Grid.Volume/dt .* (por .* obj.drhodp(Index.Start:Index.End) + rho .*dpor);
                    
                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                    J_MB_P = J_MB_P + spdiags(DiagVecs, DiagIndx, N, N);
                    
                    % 2. J_PT
                    dMupx = obj.UpWind{1,1+f}.x * ( obj.dMobdT(Index.Start:Index.End, 1) .* rho + obj.Mob(Index.Start:Index.End, 1) .* obj.drhodT(Index.Start:Index.End) );
                    dMupy = obj.UpWind{1,1+f}.y * ( obj.dMobdT(Index.Start:Index.End, 1) .* rho + obj.Mob(Index.Start:Index.End, 1) .* obj.drhodT(Index.Start:Index.End) );
                    dMupz = obj.UpWind{1,1+f}.z * ( obj.dMobdT(Index.Start:Index.End, 1) .* rho + obj.Mob(Index.Start:Index.End, 1) .* obj.drhodT(Index.Start:Index.End) );
                    
                    vecX1 = min(reshape(obj.U{1,1+f}.x(1:Nx  ,:     ,:     ), N, 1), 0) .* dMupx;
                    vecX2 = max(reshape(obj.U{1,1+f}.x(2:Nx+1,:     ,:     ), N, 1), 0) .* dMupx;
                    vecY1 = min(reshape(obj.U{1,1+f}.y(:     ,1:Ny  ,:     ), N, 1), 0) .* dMupy;
                    vecY2 = max(reshape(obj.U{1,1+f}.y(:     ,2:Ny+1,:     ), N, 1), 0) .* dMupy;
                    vecZ1 = min(reshape(obj.U{1,1+f}.z(:     ,:     ,1:Nz  ), N, 1), 0) .* dMupz;
                    vecZ2 = max(reshape(obj.U{1,1+f}.z(:     ,:     ,2:Nz+1), N, 1), 0) .* dMupz;
                    acc =  Grid.Volume/dt .* por .* obj.drhodT(Index.Start:Index.End) ;
                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                    J_MB_T = spdiags(DiagVecs,DiagIndx,N,N);
            end
            
            % Add Wells
            if f == 0 % only for reservoir
                [J_MB_P, J_MB_T, ~, ~] = obj.AddWellsToJacobian(J_MB_P, J_MB_T, Medium.State, Wells, Medium.K(:,1));
            end
        end
        function [J_EB_P , J_EB_T] = BuildMediumEnergyBalanceJacobian(obj, Medium, Wells, Grid, dt, Index, f)
            % Create local variables
            N  = Grid.N;
            P   = Medium.State.Properties('P_1').Value;
            T   = Medium.State.Properties('T').Value;
            rho = Medium.State.Properties('rho_1').Value;
            h   = Medium.State.Properties('h_1').Value;
            Rho_rock = Medium.Rho;
            Medium.ComputeDerPorosity(P);
            por = Medium.Por;
            dpor = Medium.DPor;
            
            % 1.J_EB_P
            % 1.a Pressure Block
            J_EB_P = obj.Thph{1, 1+f};
            switch class(Grid)
                case('corner_point_grid')
                    N_Face = length(obj.UpWind{1,1+f});
                    nc = Grid.N;
                    nf = length(Grid.Trans);
                    C = [ Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index , Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index ];
                    D1 = [ -double(obj.U{1,1+f}>=0)+double(obj.U{1,1+f}<0)  , double(obj.U{1,1+f}>=0)-double(obj.U{1,1+f}<0)]; D1(D1==1)=0;
                    D2 = [ -double(obj.U{1,1+f}>=0)+double(obj.U{1,1+f}<0)  , double(obj.U{1,1+f}>=0)-double(obj.U{1,1+f}<0)];
                    UpwindPermutation1 = sparse([(1:nf)'; (1:nf)'], C, D1, nf, nc)';
                    UpwindPermutation2 = sparse([(1:nf)'; (1:nf)'], C, D2, nf, nc)';
                    
                    % 1.b: compressibility part
                    Mob_drhodp_dhdp = obj.Mob(Index.Start:Index.End, 1) .* ( obj.drhodp(Index.Start:Index.End) .* h + obj.dhdp(Index.Start:Index.End) .* rho );
                    dMobUpwind = obj.UpWind{1,1+f} * Mob_drhodp_dhdp;
                    FluxDerivative = dMobUpwind .* abs(obj.U{1,1+f});
                    acc = (Grid.Volume/dt) .* ( obj.Cp.*(por.*obj.drhodp(Index.Start:Index.End)+rho.*dpor) + Medium.Cpr.*(-dpor).*Rho_rock ) .* T;
                    J_EB_P = J_EB_P + spdiags(acc,0,N,N) + (UpwindPermutation1 * spdiags(FluxDerivative,0,N_Face,N_Face) * UpwindPermutation2')';
                    
                    % 2. J_EB_T
                    J_EB_T = obj.Tk{1, 1+f};
                    Mob  = obj.Mob(Index.Start:Index.End, 1);
                    dMobdT_drhodT_dhdt = obj.dMobdT(Index.Start:Index.End) .* rho .* h + obj.drhodT(Index.Start:Index.End) .* Mob .* h  + obj.dhdT(Index.Start:Index.End) .* rho .* Mob;
                    dMobUpwind = obj.UpWind{1,1+f} * dMobdT_drhodT_dhdt;
                    FluxDerivative = dMobUpwind .* abs(obj.U{1,1+f});
                    acc = (Grid.Volume/dt) .* ( obj.Cp .* por .*( obj.drhodT(Index.Start:Index.End) .* T + rho ) + Medium.Cpr.* (1-por) .* Rho_rock );
                    J_EB_T = J_EB_T + spdiags(acc,0,N,N) + (UpwindPermutation1 * spdiags(FluxDerivative,0,N_Face,N_Face) * UpwindPermutation2')';
                    
                case('cartesian_grid')
                    Nx = Grid.Nx;
                    Ny = Grid.Ny;
                    Nz = Grid.Nz;
                    % 1.b: compressibility part
                    dMupx = obj.UpWind{1,1+f}.x * ( obj.Mob(Index.Start:Index.End, 1) .* ( obj.drhodp(Index.Start:Index.End) .* h + obj.dhdp(Index.Start:Index.End) .* rho ) );
                    dMupy = obj.UpWind{1,1+f}.y * ( obj.Mob(Index.Start:Index.End, 1) .* ( obj.drhodp(Index.Start:Index.End) .* h + obj.dhdp(Index.Start:Index.End) .* rho ) );
                    dMupz = obj.UpWind{1,1+f}.z * ( obj.Mob(Index.Start:Index.End, 1) .* ( obj.drhodp(Index.Start:Index.End) .* h + obj.dhdp(Index.Start:Index.End) .* rho ) );
                    
                    vecX1 = min(reshape(obj.U{1,1+f}.x(1:Nx  ,:     ,:     ), N, 1), 0) .* dMupx;
                    vecX2 = max(reshape(obj.U{1,1+f}.x(2:Nx+1,:     ,:     ), N, 1), 0) .* dMupx;
                    vecY1 = min(reshape(obj.U{1,1+f}.y(:     ,1:Ny  ,:     ), N, 1), 0) .* dMupy;
                    vecY2 = max(reshape(obj.U{1,1+f}.y(:     ,2:Ny+1,:     ), N, 1), 0) .* dMupy;
                    vecZ1 = min(reshape(obj.U{1,1+f}.z(:     ,:     ,1:Nz  ), N, 1), 0) .* dMupz;
                    vecZ2 = max(reshape(obj.U{1,1+f}.z(:     ,:     ,2:Nz+1), N, 1), 0) .* dMupz;
                    acc = (Grid.Volume/dt) .* ( obj.Cp.*(por.*obj.drhodp(Index.Start:Index.End)+rho.*dpor) + Medium.Cpr.*(-dpor).*Rho_rock ) .* T;
                    
                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                    J_EB_P = J_EB_P + spdiags(DiagVecs, DiagIndx, N, N);
                    
                    % 2. J_EB_T
                    J_EB_T = obj.Tk{1, 1+f};
                    Mob  = obj.Mob(Index.Start:Index.End, 1);
                    
                    dMupx = obj.UpWind{1,1+f}.x * (obj.dMobdT(Index.Start:Index.End) .* rho .* h + obj.drhodT(Index.Start:Index.End) .* Mob .* h  + obj.dhdT(Index.Start:Index.End) .* rho .* Mob);
                    dMupy = obj.UpWind{1,1+f}.y * (obj.dMobdT(Index.Start:Index.End) .* rho .* h + obj.drhodT(Index.Start:Index.End) .* Mob .* h  + obj.dhdT(Index.Start:Index.End) .* rho .* Mob);
                    dMupz = obj.UpWind{1,1+f}.z * (obj.dMobdT(Index.Start:Index.End) .* rho .* h + obj.drhodT(Index.Start:Index.End) .* Mob .* h  + obj.dhdT(Index.Start:Index.End) .* rho .* Mob);
                    % Construct J_EB_T block
                    vecX1 = min(reshape(obj.U{1,1+f}.x(1:Nx  ,:     ,:     ), N, 1), 0) .* dMupx;
                    vecX2 = max(reshape(obj.U{1,1+f}.x(2:Nx+1,:     ,:     ), N, 1), 0) .* dMupx;
                    vecY1 = min(reshape(obj.U{1,1+f}.y(:     ,1:Ny  ,:     ), N, 1), 0) .* dMupy;
                    vecY2 = max(reshape(obj.U{1,1+f}.y(:     ,2:Ny+1,:     ), N, 1), 0) .* dMupy;
                    vecZ1 = min(reshape(obj.U{1,1+f}.z(:     ,:     ,1:Nz  ), N, 1), 0) .* dMupz;
                    vecZ2 = max(reshape(obj.U{1,1+f}.z(:     ,:     ,2:Nz+1), N, 1), 0) .* dMupz;
                    acc = (Grid.Volume/dt) .* ( obj.Cp .* por .*( obj.drhodT(Index.Start:Index.End) .* T + rho ) + Medium.Cpr.* (1-por) .* Rho_rock );
                    
                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                    J_EB_T = J_EB_T + spdiags(DiagVecs, DiagIndx, N, N);
            end
            
            % Add Wells
            if f == 0 % only for reservoir
                [~, ~, J_EB_P, J_EB_T] = obj.AddWellsToJacobian(J_EB_P, J_EB_T, Medium.State, Wells, Medium.K(:,1));
            end
        end
        function JacobianFull = BuildFullJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            %% Jacobian's assembly for single phase geothermal reservoir
            % | J_MB_P   J_MB_T |  |dP|      | R_MB |
            % | J_EB_P   J_EB_T |  |dT|  =  -| R_EB |
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            Reservoir = ProductionSystem.Reservoir;
            Fractures = ProductionSystem.FracturesNetwork.Fractures;
            Wells = ProductionSystem.Wells;       
            
            J_MB_P = [];  J_MB_T = [];
            J_EB_P = [];  J_EB_T = [];
            for ph=1:obj.NofPhases
                %% Jacobian of the reservoir 
                Index.Start = 1;
                Index.End = Nm;
                % Mass Balance Jacobian blocks
                [J_MB_P_Reservoir, J_MB_T_Reservoir] = BuildMediumMassBalanceJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0); % 0 = reservoir
                % Energy Balance Jacobian blocks
                [J_EB_P_Reservoir, J_EB_T_Reservoir] = BuildMediumEnergyBalanceJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0); % 0 = reservoir
                
                %% Jacobian of the fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Nf = DiscretizationModel.FracturesGrid.N;
                    Index.Start = DiscretizationModel.Index_Local_to_Global( Nm , f , 1     );
                    Index.End   = DiscretizationModel.Index_Local_to_Global( Nm , f , Nf(f) );
                    % Mass Balance Jacobian blocks
                    [J_MB_P_Fractures, J_MB_T_Fractures] = BuildMediumMassBalanceJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f);
                    J_MB_P  = blkdiag(J_MB_P, J_MB_P_Fractures);
                    J_MB_T  = blkdiag(J_MB_T, J_MB_T_Fractures);
                    % Energy Balance Jacobian blocks
                    [J_EB_P_Fractures, J_EB_T_Fractures] = BuildMediumEnergyBalanceJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f);
                    J_EB_P  = blkdiag(J_EB_P, J_EB_P_Fractures);
                    J_EB_T  = blkdiag(J_EB_T, J_EB_T_Fractures);
                end
                
                % Adding the blocks of matrix Jacobian to the beginning of the FullJacobian
                J_MB_P  = blkdiag(J_MB_P_Reservoir, J_MB_P);
                J_MB_T  = blkdiag(J_MB_T_Reservoir, J_MB_T);
                J_EB_P  = blkdiag(J_EB_P_Reservoir, J_EB_P);
                J_EB_T  = blkdiag(J_EB_T_Reservoir, J_EB_T);

                %% ADD frac-matrix and frac-frac connections in Jacobian blocks
                % Global variables
                if ProductionSystem.FracturesNetwork.Active
                    FineGrid = [DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids];
                else
                    FineGrid = DiscretizationModel.ReservoirGrid;
                end
                P = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'P_'); % useful for cross connections assembly
                rho = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'rho_'); % useful for cross connections assembly
                h = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'h_'); % useful for cross connections assembly
                for c = 1:length(DiscretizationModel.CrossConnections)
                    if isempty(DiscretizationModel.CrossConnections(c).Cells),  continue;  end
                    T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                    T_Geo_Cond = DiscretizationModel.CrossConnections(c).T_Geo_Cond;
                    UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                    i = c + Nm;
                    j = DiscretizationModel.CrossConnections(c).Cells;
                    
                    %% J_MB_P Coupling
                    J_MB_P_1_conn = - T_Geo .*  UpWind(:, ph) .* obj.Mob(j, ph) .* ( rho(j,ph) + obj.drhodp(j, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_MB_P_2_conn = - T_Geo .* ~UpWind(:, ph) .* obj.Mob(i, ph) .* ( rho(i,ph) + obj.drhodp(i, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_MB_P_conn  = J_MB_P_1_conn + J_MB_P_2_conn;
                    
                    % frac - mat or frac1 - frac2
                    if sum(J_MB_P(i, j))~=0 || sum(J_MB_P(j, i))~=0
                        error('J_MB_P(i, j) or J_MB_P(j, i) is not zero!');
                    end
                    J_MB_P(i, j) = J_MB_P_conn;
                    J_MB_P(i, i) = J_MB_P(i, i) - sum(J_MB_P_conn);
                    % mat-frac or frac2 - frac1
                    J_MB_P(j, i) = J_MB_P_conn';
                    J_MB_P(sub2ind([Nt, Nt], j, j)) = J_MB_P(sub2ind([Nt, Nt], j, j)) - J_MB_P_conn;
                    
                    %% J_MB_T Coupling
                    J_MB_T_1_conn = - T_Geo .*  UpWind(:, ph) .* (P(j, ph) - P(i, ph)).* ( rho(j, ph) .* obj.dMobdT(j, ph) + obj.Mob(j, ph) .* obj.drhodT(j, ph) );
                    J_MB_T_2_conn = - T_Geo .* ~UpWind(:, ph) .* (P(j, ph) - P(i, ph)).* ( rho(i, ph) .* obj.dMobdT(i, ph) + obj.Mob(i, ph) .* obj.drhodT(i, ph) );
                    
                    % frac - mat or frac1 - frac2
                    if sum(J_MB_T(i, j))~=0 || sum(J_MB_T(j, i))~=0
                        error('J_MB_T(i, j) or J_MB_T(j, i) is not zero!');
                    end
                    J_MB_T(i, j) = J_MB_T_1_conn;
                    J_MB_T(i, i) = J_MB_T(i, i) + sum(J_MB_T_2_conn);
                    % mat-frac or frac2 - frac1
                    J_MB_T(j, i) = - J_MB_T_2_conn;
                    % diag of mat or frac2
                    J_MB_T(sub2ind([Nt, Nt], j, j)) = J_MB_T(sub2ind([Nt, Nt], j, j)) - J_MB_T_1_conn;
                    
                    %% J_EB_P Coupling
                    J_EB_P_1_conn = - T_Geo .*  UpWind(:, ph) .* obj.Mob(j, ph) .* ( (P(j, ph) - P(i, ph)) .* (rho(j, ph) .* obj.dhdp(j, ph) + h(j, ph) .* obj.drhodp(j, ph)) ...
                        + rho(j, ph) .* h(j, ph) );
                    J_EB_P_2_conn = - T_Geo .* ~UpWind(:, ph) .* obj.Mob(i, ph) .* ( (P(j, ph) - P(i, ph)) .* (rho(i, ph) .* obj.dhdp(i, ph) + h(i, ph) .* obj.drhodp(i, ph)) ...
                        + rho(i, ph) .* h(i, ph) );
                    J_EB_P_conn = J_EB_P_1_conn + J_EB_P_2_conn;
                    % frac - mat or frac1 - frac2
                    if sum(J_EB_P(i, j))~=0 || sum(J_EB_P(j, i))~=0
                        error('J_EB_P(i, j) or J_EB_P(j, i) is not zero!');
                    end
                    J_EB_P(i, j) = J_EB_P_conn;
                    J_EB_P(i, i) = J_EB_P(i, i) - sum(J_EB_P_conn);
                    % mat-frac or frac2 - frac1
                    J_EB_P(j, i) = J_EB_P_conn';
                    % diag of mat or frac2
                    J_EB_P(sub2ind([Nt, Nt], j, j)) = J_EB_P(sub2ind([Nt, Nt], j, j)) - J_EB_P_conn;
                    
                    %% J_EB_T Coupling
                    J_EB_T_1_conn = - T_Geo .*  UpWind(:, ph) .* ( (P(j, ph) - P(i, ph)) .* ( obj.dMobdT(j, ph) .* rho(j, ph) .* h(j, ph) + obj.Mob(j, ph) .* obj.drhodT(j, ph) .* h(j, ph) ...
                        + obj.Mob(j, ph) .* rho(j, ph) .* obj.dhdT(j, ph) ) );
                    J_EB_T_2_conn = - T_Geo .* ~UpWind(:, ph) .* ( (P(j, ph) - P(i, ph)) .* ( obj.dMobdT(i, ph) .* rho(i, ph) .* h(i, ph) + obj.Mob(i, ph) .* obj.drhodT(i, ph) .* h(i, ph) ...
                        + obj.Mob(i, ph) .* rho(i, ph) .* obj.dhdT(i, ph) ) );
                    % frac - mat or frac1 - frac2
                    if sum(J_EB_T(i, j))~=0 || sum(J_EB_T(j, i))~=0
                        error('J_EB_T(i, j) or J_EB_T(j, i) is not zero!');
                    end
                    J_EB_T(i, j) = J_EB_T_1_conn';
                    J_EB_T(i, i) = J_EB_T(i, i) + sum(J_EB_T_2_conn);
                    % mat-frac or frac2 - frac1
                    J_EB_T(j, i) = -J_EB_T_2_conn;
                    % diag of mat or frac2
                    J_EB_T(sub2ind([Nt, Nt], j, j)) = J_EB_T(sub2ind([Nt, Nt], j, j)) - J_EB_T_1_conn;
                end
            end
            
            % Build & Stack Jacobian
            JacobianFull = [J_MB_P, J_MB_T ; J_EB_P, J_EB_T];
        end
        function ConstrainedPressureResidual(obj, FluidModel, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute source terms
            [WellMassFlux, ~] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            FracturesMassFlux = zeros(DiscretizationModel.N, obj.NofPhases);
            if ProductionSystem.FracturesNetwork.Active
                [FracturesMassFlux, ~, ~] = obj.ComputeSourceTerms_frac_mat(ProductionSystem, DiscretizationModel);
            end
            % Initialise residual vector (Nph * N, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;
            if ProductionSystem.FracturesNetwork.Active
                Nf = DiscretizationModel.FracturesGrid.N;
            else
                Nf = 0;
            end
            Nt = DiscretizationModel.N;
            Residual_MB = zeros( Nt , 1 );
            
            %% Computing Transmissibilities for reservoir
            Index.Start = 1;
            Index.End = Nm;
            for ph=1:obj.NofPhases
                % Compute mass flow transmissibility and gravity
                [obj.Tph{ph, 1}, obj.Gph{ph, 1}] = obj.MatrixAssembler.TransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
            end
            
            %% Computing Transmissibilities for fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Index.Start = Index.End+1;
                Index.End = Index.Start + Nf(f) - 1;
                for ph=1:obj.NofPhases
                    % Compute mass flow transmissibility and gravity
                    [obj.Tph{ph, 1+f}, obj.Gph{ph, 1+f}] = obj.MatrixAssembler.TransmissibilityMatrix( ...
                        DiscretizationModel.FracturesGrid.Grids(f), ...
                        obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph), ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value, ...
                        obj.GravityModel.RhoInt{ph, 1+f} );
                end
            end
            
            %% Computing Mass Balance Residual
            % Reservoir
            Index.Start = 1;
            Index.End = Nm;
            Residual_MB_Reservoir = BuildMediumMassBalanceResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, WellMassFlux, FracturesMassFlux, 0);
            Residual_MB((ph-1)*Nt + Index.Start: (ph-1)*Nt + Index.End) = Residual_MB_Reservoir;
            % Fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Index.Start = Index.End+1;
                Index.End = Index.Start + Nf(f) - 1;
                Residual_MB_Fractures = BuildMediumMassBalanceResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, WellMassFlux, FracturesMassFlux, f);
                Residual_MB(Index.Start:Index.End) = Residual_MB_Fractures;
            end
            
            %% Jacobian's assembly for single phase geothermal reservoir
            % |J_MB_P| |dP| = -|R_MB|
            J_MB_P = [];
            for ph=1:obj.NofPhases
                %% Jacobian of the reservoir 
                Index.Start = 1;
                Index.End = Nm;
                % Flow Jacobian blocks
                [J_MB_P_Reservoir, ~] = BuildMediumMassBalanceJacobian(obj, ProductionSystem.Reservoir, ProductionSystem.Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0); % 0 = reservoir
                
                %% Jacobian of the fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Nf = DiscretizationModel.FracturesGrid.N;
                    Index.Start = DiscretizationModel.Index_Local_to_Global( Nm , f , 1     );
                    Index.End   = DiscretizationModel.Index_Local_to_Global( Nm , f , Nf(f) );
                    % Flow Jacobian Block
                    [J_MB_P_Fractures, ~] = BuildMediumMassBalanceJacobian(obj, ProductionSystem.FracturesNetwork.Fractures(f), ProductionSystem.Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f);
                    J_MB_P  = blkdiag(J_MB_P, J_MB_P_Fractures);
                end
                % Adding the blocks of matrix Jacobian to the beginning of the FullJacobian
                J_MB_P  = blkdiag(J_MB_P_Reservoir, J_MB_P);
                
                %% ADD frac-matrix and frac-frac connections in Jacobian blocks
                % Global variables
                if ProductionSystem.FracturesNetwork.Active
                    FineGrid = [DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids];
                else
                    FineGrid = DiscretizationModel.ReservoirGrid;
                end
                P = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'P_'); % useful for cross connections assembly
                rho = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'rho_'); % useful for cross connections assembly
                for c = 1:length(DiscretizationModel.CrossConnections)
                    T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                    UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                    i = c + Nm;
                    j = DiscretizationModel.CrossConnections(c).Cells;
                    
                    % J_PP Coupling
                    J_MB_P_1_conn = - T_Geo .*  UpWind(:, ph) .* obj.Mob(j, ph) .* ( rho(j,ph) + obj.drhodp(j, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_MB_P_2_conn = - T_Geo .* ~UpWind(:, ph) .* obj.Mob(i, ph) .* ( rho(i,ph) + obj.drhodp(i, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_MB_P_conn  = J_MB_P_1_conn + J_MB_P_2_conn;
                    % frac - mat or frac1 - frac2
                    if sum(J_MB_P(i, j))~=0 || sum(J_MB_P(j, i))~=0
                        error('J_MB_P(i, j) or J_MB_P(j, i) is not zero!');
                    end
                    J_MB_P(i, j) = J_MB_P_conn;
                    J_MB_P(i, i) = J_MB_P(i, i) - sum(J_MB_P_conn);
                    % mat-frac or frac2 - frac1
                    J_MB_P(j, i) = J_MB_P_conn';
                    J_MB_P(sub2ind([Nt, Nt], j, j)) = J_MB_P(sub2ind([Nt, Nt], j, j)) - J_MB_P_conn;
                end
            end
            
            %% Solve for deltaP
            if isprop(DiscretizationModel,'OperatorsHandler')
                deltaP_Coarse = (  DiscretizationModel.OperatorsHandler.ADMRest * J_MB_P * DiscretizationModel.OperatorsHandler.ADMProl{1} ) \ ...
                                (- DiscretizationModel.OperatorsHandler.ADMRest * Residual_MB );
                deltaP = DiscretizationModel.OperatorsHandler.ADMProl{1} * deltaP_Coarse;
            else
                deltaP = J_MB_P\(-Residual_MB);
            end
            
            %% Update Reservoir State
            Pm = ProductionSystem.Reservoir.State.Properties('P_1');
            Pm.update(deltaP(1:Nm));

            %% Update fracture state
            if ProductionSystem.FracturesNetwork.Active
                EP = Nm;
                Nf = DiscretizationModel.FracturesGrid.N;
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    IP = EP+1;
                    EP = IP + Nf(f) - 1;
                    Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_1']);
                    Pf.update(deltaP(IP:EP));
                end
            end
            
            %% Update the peoperties and well state
            obj.ComputeProperties(ProductionSystem, FluidModel);
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
        end
        function ConstrainedTemperatureResidual(obj, FluidModel, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute source terms
            [~, WellHeatConvectionFlux] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            FracturesHeatConvectionFlux = zeros(DiscretizationModel.N, obj.NofPhases); % Heat Convection flux betweem each two media
            FracturesHeatConductionFlux = zeros(DiscretizationModel.N, 1);             % Heat Conduction flux betweem each two media
            if ProductionSystem.FracturesNetwork.Active
                [~, FracturesHeatConvectionFlux, FracturesHeatConductionFlux] = obj.ComputeSourceTerms_frac_mat(ProductionSystem, DiscretizationModel);
            end
            FracturesHeatConductionFlux = zeros(DiscretizationModel.N, 1);
            
            % Initialise residual vector (Nph * N, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;
            if ProductionSystem.FracturesNetwork.Active
                Nf = DiscretizationModel.FracturesGrid.N;
            else
                Nf = 0;
            end
            Nt = DiscretizationModel.N;
            Residual_EB = zeros( Nt , 1 );
            
            %% Computing Transmissibilities for reservoir
            Index.Start = 1;
            Index.End = Nm;
            for ph=1:obj.NofPhases
                % Compute heat convection transmissibility and gravity
                [obj.Thph{ph, 1}, obj.Ghph{ph, 1}] = obj.MatrixAssembler.ConvectiveHeatTransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    ProductionSystem.Reservoir.State.Properties(['h_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
            end
            % Compute heat conduction transmissibility
            obj.Tk{1,1} = obj.MatrixAssembler.ConductiveHeatTransmissibilityMatrix( DiscretizationModel.ReservoirGrid , ProductionSystem.Reservoir.State );
            
            %% Computing Transmissibilities for fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Index.Start = Index.End+1;
                Index.End = Index.Start + Nf(f) - 1;
                for ph=1:obj.NofPhases
                    % Compute heat convection transmissibility and gravity
                    [obj.Thph{ph, 1+f}, obj.Ghph{ph, 1+f}] = obj.MatrixAssembler.ConvectiveHeatTransmissibilityMatrix( ...
                        DiscretizationModel.FracturesGrid.Grids(f), ...
                        obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph), ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value, ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['h_',num2str(ph)]).Value, ...
                        obj.GravityModel.RhoInt{ph, 1+f} );
                end
                % Compute heat conduction transmissibility
                obj.Tk{1,1+f} = obj.MatrixAssembler.ConductiveHeatTransmissibilityMatrix( DiscretizationModel.FracturesGrid.Grids(f) , ProductionSystem.FracturesNetwork.Fractures(f).State );
            end

            for ph=1:obj.NofPhases
                %% Energy Balance Residual
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                Residual_EB_Reservoir = BuildMediumHeatResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, WellHeatConvectionFlux, FracturesHeatConvectionFlux, FracturesHeatConductionFlux, 0, ph);
                Residual_EB(Index.Start: Index.End) = Residual_EB_Reservoir;
                % Fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End + 1;
                    Index.End = Index.Start + Nf(f) - 1;
                    Residual_EB_Fractures = BuildMediumHeatResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, WellHeatConvectionFlux, FracturesHeatConvectionFlux, FracturesHeatConductionFlux, f, ph);
                    Residual_EB(Index.Start:Index.End) = Residual_EB_Fractures; 
                end
            end
            
            %% Jacobian's assembly for single phase geothermal reservoir
            % |J_EB_T| |dT| = -|R_EB|
            J_EB_T = [];
            for ph=1:obj.NofPhases
                %% Jacobian of the reservoir 
                Index.Start = 1;
                Index.End = Nm;
                % Heat Jacobian blocks
                [~, J_EB_T_Reservoir] = BuildMediumHeatJacobian(obj, ProductionSystem.Reservoir, ProductionSystem.Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0, ph); % 0 = reservoir
                
                %% Jacobian of the fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Nf = DiscretizationModel.FracturesGrid.N;
                    Index.Start = DiscretizationModel.Index_Local_to_Global( Nm , f , 1     );
                    Index.End   = DiscretizationModel.Index_Local_to_Global( Nm , f , Nf(f) );
                    % Heat Jacobian Block
                    [~, J_EB_T_Fractures] = BuildMediumHeatJacobian(obj, ProductionSystem.FracturesNetwork.Fractures(f), ProductionSystem.Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f, ph);
                    J_EB_T  = blkdiag(J_EB_T, J_EB_T_Fractures);
                end
                % Adding the blocks of matrix Jacobian to the beginning of the FullJacobian
                J_EB_T  = blkdiag(J_EB_T_Reservoir, J_EB_T);
                
                %% ADD frac-matrix and frac-frac connections in Jacobian blocks
                % Global variables
                if ProductionSystem.FracturesNetwork.Active
                    FineGrid = [DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids];
                else
                    FineGrid = DiscretizationModel.ReservoirGrid;
                end
                P   = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'P_'); % useful for cross connections assembly
                rho = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'rho_'); % useful for cross connections assembly
                h   = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'h_'); % useful for cross connections assembly
                for c = 1:length(DiscretizationModel.CrossConnections)
                    T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                    T_Geo_Cond = DiscretizationModel.CrossConnections(c).T_Geo_Cond;
                    UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                    i = c + Nm;
                    j = DiscretizationModel.CrossConnections(c).Cells;
                    
                    %% J_EB_T Coupling
                    J_EB_T_1_conn = - T_Geo .*  UpWind(:, ph) .* ( (P(j, ph) - P(i, ph)) .* ( obj.dMobdT(j, ph) .* rho(j, ph) .* h(j, ph) + obj.Mob(j, ph) .* obj.drhodT(j, ph) .* h(j, ph) ...
                        + obj.Mob(j, ph) .* rho(j, ph) .* obj.dhdT(j, ph) ) ); - T_Geo_Cond .*  UpWind(:, ph);
                    J_EB_T_2_conn = - T_Geo .* ~UpWind(:, ph) .* ( (P(j, ph) - P(i, ph)) .* ( obj.dMobdT(i, ph) .* rho(i, ph) .* h(i, ph) + obj.Mob(i, ph) .* obj.drhodT(i, ph) .* h(i, ph) ...
                        + obj.Mob(i, ph) .* rho(i, ph) .* obj.dhdT(i, ph) ) ); - T_Geo_Cond .*  UpWind(:, ph);
                    % frac - mat or frac1 - frac2
                    if sum(J_EB_T(i, j))~=0 || sum(J_EB_T(j, i))~=0
                        error('J_EB_T(i, j) or J_EB_T(j, i) is not zero!');
                    end
                    J_EB_T(i, j) = J_EB_T_1_conn';
                    J_EB_T(i, i) = J_EB_T(i, i) + sum(J_EB_T_2_conn);
                    % mat-frac or frac2 - frac1
                    J_EB_T(j, i) = -J_EB_T_2_conn;
                    % diag of mat or frac2
                    J_EB_T(sub2ind([Nt, Nt], j, j)) = J_EB_T(sub2ind([Nt, Nt], j, j)) - J_EB_T_1_conn;
                end
            end
            
            %% Solve for deltaT
            if isprop(DiscretizationModel,'OperatorsHandler')
                deltaTc = (  DiscretizationModel.OperatorsHandler.ADMRest * J_EB_T * DiscretizationModel.OperatorsHandler.ADMProl{2} ) \ ...
                          (- DiscretizationModel.OperatorsHandler.ADMRest * Residual_EB );
                deltaT = DiscretizationModel.OperatorsHandler.ADMProl{2} * deltaTc;
            else
                deltaT = J_EB_T\(-Residual_EB);
            end
            
            %% Update Reservoir State
            Tm = ProductionSystem.Reservoir.State.Properties('T');
            Tm.update(deltaT(1:Nm));

            %% Update fracture state
            % Update Pressure
            if ProductionSystem.FracturesNetwork.Active
                EP = Nm;
                Nf = DiscretizationModel.FracturesGrid.N;
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    IP = EP+1;
                    EP = IP + Nf(f) - 1;
                    Tf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('T');
                    Tf.update(deltaT(IP:EP));
                end
            end
            
            %% Update the peoperties and well state
            obj.ComputeProperties(ProductionSystem, FluidModel);
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
        end
        function delta = UpdateState(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                return
            else
                % Inflexion point correction
%                 delta = obj.UpdateTemperature(delta, ProductionSystem, FluidModel, DiscretizationModel);
                
                Nm = DiscretizationModel.ReservoirGrid.N;
                Nt = DiscretizationModel.N;
                deltaP = delta(1:Nt);
                deltaT = delta(Nt+1:2*Nt);
                
                %% Update reservoir state
                % 1. Update Pressure
                Pm = ProductionSystem.Reservoir.State.Properties('P_1');
                Pm.update(deltaP(1:Nm));
                % 2. Update Temperature
                Tm = ProductionSystem.Reservoir.State.Properties('T');
                Tm.update(deltaT(1:Nm));
            
                %% Update fracture state
                % 1. Update Pressure
                if ProductionSystem.FracturesNetwork.Active
                    EP = Nm;
                    Nf = DiscretizationModel.FracturesGrid.N;
                    for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                        IP = EP+1;
                        EP = IP + Nf(f) - 1;
                        Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('P_1');
                        Pf.update(deltaP(IP:EP));
                        % 2. Update Temperature
                        T = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('T');
                        T.update(deltaT(IP:EP));
                    end
                end
                
                % Update properties and derivatives
                obj.ComputeProperties(ProductionSystem, FluidModel)            
                % Update wells
                ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
            end
        end
        function [WellMassFlux, WellHeatConvectionFlux]= ComputeSourceTerms(obj, N, Wells)
            WellMassFlux           = zeros(N, obj.NofPhases);
            WellHeatConvectionFlux = zeros(N, obj.NofPhases);  
            % Injectors
            for i=1:Wells.NofInj
                c = Wells.Inj(i).Cells;
                WellMassFlux(c, :)           = Wells.Inj(i).QPhases(:,:);
                WellHeatConvectionFlux(c, :) = Wells.Inj(i).Qh(:,:);
            end
            % Producers
            for i=1:Wells.NofProd
                c = Wells.Prod(i).Cells;
                WellMassFlux(c, :)           = Wells.Prod(i).QPhases(:,:);
                WellHeatConvectionFlux(c, :) = Wells.Prod(i).Qh(:,:);
            end
        end
        function [FracturesMassFlux, FracturesHeatConvectionFlux, FracturesHeatConductionFlux]= ComputeSourceTerms_frac_mat(obj, ProductionSystem, DiscretizationModel)
            FracturesMassFlux           = zeros(DiscretizationModel.N, obj.NofPhases);  % Mass Flow flux between each two media 
            FracturesHeatConvectionFlux = zeros(DiscretizationModel.N, obj.NofPhases);  % Heat Convection flux betweem each two media
            FracturesHeatConductionFlux = zeros(DiscretizationModel.N, 1);              % Heat Conduction flux betweem each two media
            Nm = DiscretizationModel.ReservoirGrid.N;

            % Global variables
            P   = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'P_'); % useful for cross connections assembly
            rho = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'rho_'); % useful for cross connections assembly
            h   = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'h_'); % useful for cross connections assembly
            T   = ProductionSystem.CreateGlobalSinglePhaseVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], 'T'); % useful for cross connections assembly
            for c = 1:length(DiscretizationModel.CrossConnections)
                j = DiscretizationModel.CrossConnections(c).Cells;
                i = c + Nm;
                T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                T_Geo_Cond = DiscretizationModel.CrossConnections(c).T_Geo_Cond;
                UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                
                for ph=1:obj.NofPhases
                    FracturesMassFlux(i, ph) = FracturesMassFlux(i, ph) + sum(T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* (P(j, ph) - P(i, ph)) + ...
                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* (P(j, ph) - P(i, ph))) );
                    FracturesMassFlux(j, ph) = FracturesMassFlux(j, ph) + T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* (P(i, ph) - P(j, ph)) + ...
                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* (P(i, ph) - P(j, ph)) );
                    
                    FracturesHeatConvectionFlux(i, ph) = FracturesHeatConvectionFlux(i, ph) + sum(T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* h(j,ph) .* (P(j, ph) - P(i, ph)) + ...
                        ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* h(i,ph) .* (P(j, ph) - P(i, ph))) );
                    FracturesHeatConvectionFlux(j, ph) = FracturesHeatConvectionFlux(j, ph) + T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* h(j,ph) .* (P(i, ph) - P(j, ph)) + ...
                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* h(i,ph) .* (P(i, ph) - P(j, ph)) );
                end
                
                FracturesHeatConductionFlux(i) = FracturesHeatConductionFlux(i) + sum( T_Geo_Cond .* (T(j) - T(i)) );
                FracturesHeatConductionFlux(j) = FracturesHeatConductionFlux(j) + T_Geo_Cond .* (T(i) - T(j));
            end
        end
        function [J_MB_P, J_MB_T] = AddWellsToMassBalanceJacobian(obj, J_MB_P, J_MB_T, State, Wells, K)
            %Injectors
            for i=1:length(Inj)
                a = Wells.Inj(i).Cells;
                dQdp  = Wells.Inj(i).ComputeWellMassFluxDerivativeWithRespectToPressure(K, obj.NofPhases);
                dQdT  = Wells.Inj(i).ComputeWellMassFluxDerivativeWithRespectToTemperature(obj.NofPhases);
                for ph=1:obj.NofPhases
                    for j=1:length(a)
                        % Add derivative of injection well mass flux with respect to pressure, to J_MB_P
                        J_MB_P(a(j),a(j)) = J_MB_P(a(j),a(j)) - dQdp(j, ph);
                        % Add derivative of injection well mass flux with respect to temperature, to J_MB_T
                        J_MB_T(a(j),a(j)) = J_MB_T(a(j),a(j)) - dQdT(j, ph);
                    end
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Wells.Prod(i).Cells;
                dQdp = Wells.Prod(i).ComputeWellMassFluxDerivativeWithRespectToPressure   (State, K, obj.Mob, obj.dMobdp, obj.drhodp, obj.NofPhases);
                dQdT = Wells.Prod(i).ComputeWellMassFluxDerivativeWithRespectToTemperature(State, K, obj.Mob, obj.dMobdT, obj.drhodT, obj.NofPhases);
                for ph=1:obj.NofPhases
                    for j=1:length(b)
                        % add derivative of production well mass flux with respect to pressure to J_MB_P
                        J_MB_P(b(j),b(j)) = J_MB_P(b(j),b(j)) - dQdp(j, ph);
                        % Add derivative of production well mass flux with respect to temperature, to J_MB_T
                        J_MB_T(b(j),b(j)) = J_MB_T(b(j),b(j)) - dQdT(j, ph);
                    end
                end
            end
        end
        function [J_EB_P, J_EB_H] = AddWellsToEnergyBalanceJacobian(obj, J_EB_P, J_EB_H, State, Wells, K)
            %Injectors
            for i=1:length(Inj)
                a = Wells.Inj(i).Cells;
                dQhdp = Wells.Inj(i).ComputeWellHeatFluxDerivativeWithRespectToPressure(K, obj.NofPhases);
                dQhdT = Wells.Inj(i).ComputeWellHeatFluxDerivativeWithRespectToTemperature(obj.NofPhases);
                for j=1:length(a)
                    % add derivative of inj well to Jpp
                    J_EB_P(a(j),a(j)) = J_MB_P(a(j),a(j)) - dQhdp(j, ph);
                    % add derivative of inj well to JTp
                    J_EB_P(a(j),a(j)) = AP2(a(j),a(j)) - dQhdp(j, ph);
                end
            end
        end
        function [J_MB_P, AT1, AP2, AT2] = AddWellsToJacobian(obj, AP, AT, State, Wells, K, ph)
            J_MB_P = AP; % add wells to pressure-pressure jacobian (mass balance)
            AT1 = AT; % add wells to pressure-heat jacobian (mass balance)
            AP2 = AP; % add wells to heat-pressure jacobian (energy balance)
            AT2 = AT; % add wells to heat-heat jacobian (energy balance)
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                dQdp = Inj(i).ComputeWellMassFluxDerivativeWithRespectToPressure(K, obj.NofPhases);
                dQhdp = Inj(i).ComputeWellHeatFluxDerivativeWithRespectToPressure(State, K, obj.NofPhases);
                dQdT = Inj(i).ComputeWellMassFluxDerivativeWithRespectToTemperature(obj.NofPhases);
                dQhdT = Inj(i).ComputeWellHeatFluxDerivativeWithRespectToTemperature(obj.NofPhases);
                for j=1:length(a)
                    % add derivative of inj well to Jpp
                    J_MB_P(a(j),a(j)) = J_MB_P(a(j),a(j)) - dQdp(j, ph);
                    % add derivative of inj well to JTp
                    AP2(a(j),a(j)) = AP2(a(j),a(j)) - dQhdp(j, ph);
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;                
                dQdp = Prod(i).ComputeWellMassFluxDerivativeWithRespectToPressure(State, K, obj.Mob, obj.drhodp, obj.dMobdp, obj.NofPhases); 
                dQdT = Prod(i).ComputeWellMassFluxDerivativeWithRespectToTemperature(State, K, obj.Mob, obj.dMobdT, obj.drhodT, obj.NofPhases);
                dQhdp = Prod(i).ComputeWellHeatFluxDerivativeWithRespectToPressure(State, K, obj.Mob, obj.drhodp, obj.dMobdp, obj.NofPhases); 
                dQhdT = Prod(i).ComputeWellHeatFluxDerivativeWithRespectToTemperature(State, K, obj.Mob, obj.dMobdT, obj.drhodT, obj.dhdT, obj.NofPhases);
                for j=1:length(b)
                    % add derivative of prod well to Jpp & JpT
                    J_MB_P(b(j),b(j)) = J_MB_P(b(j),b(j)) - dQdp(j, ph);                    
                    AT1(b(j),b(j)) = AT1(b(j),b(j)) - dQdT(j, ph);
                    % add derivative of prod well to JTp & JTT
                    AP2(b(j),b(j)) = AP2(b(j),b(j)) - dQhdp(j, ph);
                    AT2(b(j),b(j)) = AT2(b(j),b(j)) - dQhdT(j, ph);
                end
            end
        end
        function d2FdT2 = ComputeSecondDerivativeHeatConvectionFlux(obj, p, T, FluidModel)
            % Viscosity
            [mu,dmudT,d2mudT2] = FluidModel.Phases.ComputeViscosity(T);
            % Mobility
            Mob = 1./mu;
            dMobdT = -dmudT./(mu.^2);
            d2MobdT2 = (2.*dmudT.^2./mu.^3) - d2mudT2./mu.^2;
            % Density
            rho = FluidModel.Phases.ComputeDensity(p,T);
            [drhodT,d2rhodT2] = FluidModel.Phases.ComputeDrhoDT(p,T);           
            % Enthalpy
            h = FluidModel.Phases.ComputeEnthalpy(p,T);
            [dhdT,d2hdT2] = FluidModel.Phases.ComputeDhDT(p,T);
            
            % F = rho.*h.*Mob;
            % dFdT = rho.*h.*obj.dMobdT + rho.*obj.dhdT.*obj.Mob + obj.drhodT.*h.*obj.Mob;
            d2FdT2 = rho.*h.*d2MobdT2 + rho.*dhdT.*dMobdT + drhodT.*h.*dMobdT + ...
                     rho.*dhdT.*dMobdT + rho.*d2hdT2.*Mob + drhodT.*dhdT.*Mob + ...
                     drhodT.*h.*dMobdT + drhodT.*dhdT.*Mob + d2rhodT2.*h.*Mob;
        end
        function delta = UpdateTemperature(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            Nt = DiscretizationModel.N;
            deltaP_Old = delta(   1:  Nt);
            deltaT_Old = delta(Nt+1:2*Nt);
            if ProductionSystem.FracturesNetwork.Active
                    FineGrid = [DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids];
                else
                    FineGrid = DiscretizationModel.ReservoirGrid;
            end
            
            p     = ProductionSystem.CreateGlobalVariables(FineGrid, 1, 'P_');
            T_old = ProductionSystem.CreateGlobalSinglePhaseVariables(FineGrid, 'T');
            T_new = T_old + deltaT_Old;
            
            % Flux Correction - Patrick
            d2FdT2_old = obj.ComputeSecondDerivativeHeatConvectionFlux(p, T_old, FluidModel);
            d2FdT2_new = obj.ComputeSecondDerivativeHeatConvectionFlux(p, T_new, FluidModel);
            T_new = T_new.*(d2FdT2_new.*d2FdT2_old >= 0) + 0.5*(T_new + T_old).*(d2FdT2_new.*d2FdT2_old < 0);
            deltaT_New = T_new - T_old;
            
            deltaP_New = deltaP_Old .* (deltaT_New./deltaT_Old);
            delta(   1:  Nt) = deltaP_New;
            delta(Nt+1:2*Nt) = deltaT_New;
        end
        function ComputeTotalFluxes(obj, ProductionSystem, DiscretizationModel)
            % this is virtual call
        end
        function AverageMassOnCoarseBlocks(obj, ProductionSystem, FineGrid, FluidModel, ADMRest)
            T = ProductionSystem.CreateGlobalSinglePhaseVariables(FineGrid, 'T'); % No phase index for this variable
            rho = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'rho_');
            % Perform Average for ADM
            T_rest = ADMRest * T;
            rho_rest = ADMRest * rho;     
            T_rho_rest = ADMRest * (T .* rho);  
            Tav = ADMRest' * ( T_rest ./ sum(ADMRest, 2));
            delta_T = Tav - T;
                     
            rho_av = ADMRest' * ( rho_rest ./ sum(ADMRest, 2));
            delta_rho = rho_av - rho;

            % Updating the variables in reservoir
            Start=1;
            End = FineGrid(1).N;
            T = ProductionSystem.Reservoir.State.Properties('T');
            rho = ProductionSystem.Reservoir.State.Properties('rho_1');
            T.update(delta_T(Start:End));
            rho.update(delta_rho(Start:End));
            
            % Updating the variables in fractures
            for frac = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                Start = End + 1;
                End = Start + FineGrid(frac+1).N - 1;
                T = ProductionSystem.FracturesNetwork.Fractures(frac).State.Properties('T');
                rho = ProductionSystem.FracturesNetwork.Fractures(frac).State.Properties('rho_1');
                T.update(delta_T(Start:End));
                rho.update(delta_rho(Start:End));
            end
        end
        function CFL = ComputeCFLNumber(obj, ProductionSystem, DiscretizationModel, dt)
            CFL = 0;
        end
    end
end