% Geothermal MultiPhase formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Geothermal_MultiPhase_formulation < formulation
    properties
        drhodT
        d2rhodT2
        dmudT
        d2mudT2
        dhdp
        dhdT
        d2hdT2
        dMobdT
        d2MobdT2
        Kf % fluid thermal conductivity
        Tk % transmisibility of rock conductivity
        Th % transmisibility of rho .* h
    end
    methods
        function obj = Geothermal_MultiPhase_formulation(NofPhases)
            obj@formulation();
            obj.Tph = cell(NofPhases,1);
            obj.Gph = cell(NofPhases,1);
        end
        function x = GetPrimaryUnknowns(obj, ProductionSystem, DiscretizationModel)
             Nt = DiscretizationModel.N; % total grid
             Nm = DiscretizationModel.ReservoirGrid.N; % matrix grid
             x = zeros(2 * Nt, 1); % 2 unknowns: P, H
             Index.Start = 1;
             Index.End = Nm;
                     
             for ph=1:obj.NofPhases
                 % Get matrix unknowns
                 x(Index.Start     :Index.End       ) = ProductionSystem.Reservoir.State.Properties(['P_',num2str(ph)]).Value;
                 x(Nt+Index.Start  :Nt + Index.End  ) = ProductionSystem.Reservoir.State.Properties('H').Value;
                 
                 % Get fracture unknowns
                 for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                     Index.Start = Index.End+1;
                     Index.End = Index.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                     x(Index.Start     :Index.End       ) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_',num2str(ph)]).Value;
                     x(Nt+Index.Start  :Nt + Index.End  ) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('H').Value;   
                 end
             end
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            %% 1. Geothermal Properties
            FluidModel.GetTableIndex(ProductionSystem.Reservoir.State);
            
            FluidModel.GetPhaseDensities(ProductionSystem.Reservoir.State);
            FluidModel.GetPhaseSaturations(ProductionSystem.Reservoir.State);
            FluidModel.GetPhaseViscosities(ProductionSystem.Reservoir.State);
            FluidModel.GetPhaseInternalEnergies(ProductionSystem.Reservoir.State);
            FluidModel.GetTemperature(ProductionSystem.Reservoir.State);
            FluidModel.GetPhaseThermalConductivities(ProductionSystem.Reservoir.State);

            %% 2. Reservoir Properties and Derivatives
            obj.dconddp = FluidModel.ComputeDcondDp(ProductionSystem.Reservoir.State);
% you need to put obj. in front of the derivatives as they are stored as
% properties which are called when building the Jacobian matrix (similar as
% Pindex, Hindex in fluidmodel)!
            obj.drhodp = FluidModel.ComputeDrhoDp();
            [obj.drhodT,obj.d2rhodT2] = FluidModel.ComputeDrhoDT(ProductionSystem.Reservoir.State);
            [obj.dmudT,obj.d2mudT2] = FluidModel.ComputeDmuDT(ProductionSystem.Reservoir.State);
            obj.dhdp = FluidModel.ComputeDhDp(ProductionSystem.Reservoir.State);
            [obj.dhdT,obj.d2hdT2] = FluidModel.ComputeDhDT(ProductionSystem.Reservoir.State);
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.Properties('mu_1').Value); 
            [obj.dMobdT,obj.d2MobdT2] = FluidModel.ComputeDMobdT(ProductionSystem.Reservoir.State);
            %% 3. Fractures Properties and Derivatives
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
                obj.Mob = [obj.Mob; FluidModel.ComputePhaseMobilities(ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('mu_1').Value)];
                [dMobdT_f,d2MobdT2_f] = FluidModel.ComputeDMobdT(ProductionSystem.FracturesNetwork.Fractures(f).State);
                obj.dMobdT = [obj.dMobdT; dMobdT_f];
                obj.d2MobdT2 = [obj.d2MobdT2; d2MobdT2_f];
            end
        end
        %% Methods for FIM Coupling
        function Residual_P   = BuildMediumFlowResidual(obj, Medium, Grid, dt, State0, Index, qw, qf, f, ph)
            % Create local variables
            rho_old = State0.Properties(['rho_', num2str(ph)]).Value(Index.Start:Index.End);
            P_old = State0.Properties(['P_', num2str(ph)]).Value(Index.Start:Index.End); % P at previous time step
            P_new = Medium.State.Properties(['P_', num2str(ph)]).Value;
            rho_new = Medium.State.Properties(['rho_', num2str(ph)]).Value;
            depth = Grid.Depth;
            Medium.ComputePorosity(P_old);
            pv_old = Medium.Por*Grid.Volume;
            Medium.ComputePorosity(P_new);
            pv_new = Medium.Por*Grid.Volume;
            
            % Accumulation Term
            Accumulation = (pv_new.*rho_new - pv_old.*rho_old)/dt;

            % RESIDUAL
            Residual_P  = Accumulation ...
                + obj.Tph{ph, 1+f} * P_new...
                - obj.Gph{ph, 1+f} * depth...
                - qw(Index.Start:Index.End, ph)...
                - qf(Index.Start:Index.End, ph);
        end
        function Residual_T   = BuildMediumHeatResidual(obj, Medium, Grid, dt, State0, Index, qhw, qhf, RTf, f, ph)
            % Create local variables
            rho_old = State0.Properties(['rho_', num2str(ph)]).Value(Index.Start:Index.End);
            P_old = State0.Properties(['P_', num2str(ph)]).Value(Index.Start:Index.End); % P at previous time step
            P_new = Medium.State.Properties(['P_', num2str(ph)]).Value;
            T_old = State0.Properties('T').Value(Index.Start:Index.End); % T at previous time step
            T_new = Medium.State.Properties('T').Value;
            rho_new = Medium.State.Properties(['rho_', num2str(ph)]).Value;
            Rho_rock = Medium.Rho;
            depth = Grid.Depth;
            
            % Pore Volume & Rock Volume
            Medium.ComputePorosity(P_old);
            pv_old = Medium.Por*Grid.Volume;         % Old pore Volume
            mv_old = (1 - Medium.Por) * Grid.Volume; % Old rock volume
            Medium.ComputePorosity(P_new);           % Updating porosity
            pv_new = Medium.Por*Grid.Volume;         % New pore volume
            mv_new = (1 - Medium.Por) * Grid.Volume; % New rock volume
            
            % Accumulation Term
            U_eff_old = ( rho_old .* obj.Cp .* pv_old + Rho_rock .* Medium.Cpr .* mv_old ) .* T_old;
            U_eff_new = ( rho_new .* obj.Cp .* pv_new + Rho_rock .* Medium.Cpr .* mv_new ) .* T_new;
            Accumulation = (U_eff_new - U_eff_old)/dt;
            
            % RESIDUAL
            Residual_T  = Accumulation ...
                + obj.Th{ph, 1+f} * P_new ...
                - obj.Gph{ph, 1+f} * depth ...
                + obj.Tk{ph, 1+f} * T_new ...
                - qhw(Index.Start:Index.End, ph)...
                - qhf(Index.Start:Index.End, ph)...
                - RTf(Index.Start:Index.End, ph);
        end
        function ResidualFull = BuildResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute source terms
            [Qw, Qhw] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            Qf = zeros(DiscretizationModel.N, obj.NofPhases);              % Mass Flow flux between each two media 
            Qhf= zeros(DiscretizationModel.N, obj.NofPhases);              % Heat Convection flux betweem each two media
            RTf= zeros(DiscretizationModel.N, 1);                          % Heat Conduction flux betweem each two media
            if ProductionSystem.FracturesNetwork.Active
                [Qf, Qhf, RTf] = obj.ComputeSourceTerms_frac_mat(ProductionSystem, DiscretizationModel);
            end
            RTf= zeros(DiscretizationModel.N, 1);
            
            % Initialise residual vector (Nph * N, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;
            if ProductionSystem.FracturesNetwork.Active
                Nf = DiscretizationModel.FracturesGrid.N;
            else
                Nf = 0;
            end
            Nt = DiscretizationModel.N;
            ResidualFull = zeros( 2*Nt , 1 );
                        
            for ph=1:obj.NofPhases
                %% Transmissibilities of flow and heat
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                [obj.Tph{ph, 1}, obj.Gph{ph, 1}, obj.Th{ph, 1}, obj.Tk{1, 1}] = obj.TransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    ProductionSystem.Reservoir.State.Properties(['h_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
                % Fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End+1;
                    Index.End = Index.Start + Nf(f) - 1;
                    [obj.Tph{ph, 1+f}, obj.Gph{ph, 1+f}, obj.Th{ph, 1+f}, obj.Tk{1, 1+f}] = obj.TransmissibilityMatrix( ...
                        DiscretizationModel.FracturesGrid.Grids(f), ...
                        obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph), ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value, ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['h_',num2str(ph)]).Value, ...
                        obj.GravityModel.RhoInt{ph, 1+f});
                end
                
                %% Flow Residual
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                Residualm = BuildMediumFlowResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, Qw, Qf, 0, ph);
                ResidualFull((ph-1)*Nt + Index.Start: (ph-1)*Nt + Index.End) = Residualm;
                % Fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End+1;
                    Index.End = Index.Start + Nf(f) - 1;
                    Residualf = BuildMediumFlowResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, Qw, Qf, f, ph);
                    ResidualFull(Index.Start:Index.End) = Residualf;
                end
                
                %% Heat Residual
                % Reservoir
                Index.Start = Index.End + 1;
                Index.End = Index.Start + Nm - 1;
                Index_r.Start = 1;
                Index_r.End = Nm;
                Residualm = BuildMediumHeatResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index_r, Qhw, Qhf, RTf, 0, ph);
                ResidualFull(Index.Start: Index.End) = Residualm;
                % Fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End + 1;
                    Index.End = Index.Start + Nf(f) - 1;
                    Index_r.Start = Index_r.End + 1;
                    Index_r.End = Index_r.Start + Nf(f) - 1;
                    Residualf = BuildMediumHeatResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index_r, Qhw, Qhf, RTf, f, ph);
                    ResidualFull(Index.Start:Index.End) = Residualf; 
                end
            end
        end
        function [J_PP , J_PT] = BuildMediumFlowJacobian(obj, Medium, Wells, Grid, dt, Index, f, ph)
            % Create local variables
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            
            P = Medium.State.Properties(['P_', num2str(ph)]).Value;
            rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
            Medium.ComputeDerPorosity(P);
            por = Medium.Por;
            dpor = Medium.DPor;
            
            %% 1.J_PP
            % 1.a Pressure Block
            J_PP = obj.Tph{ph,1+f};
            
            % 1.b: compressibility part
            dMupx = obj.UpWind{ph,1+f}.x * ( obj.Mob(Index.Start:Index.End, ph) .* obj.drhodp(Index.Start:Index.End) );
            dMupy = obj.UpWind{ph,1+f}.y * ( obj.Mob(Index.Start:Index.End, ph) .* obj.drhodp(Index.Start:Index.End) );
            dMupz = obj.UpWind{ph,1+f}.z * ( obj.Mob(Index.Start:Index.End, ph) .* obj.drhodp(Index.Start:Index.End) );
            
            vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx  ,:     ,:     ), N, 1), 0) .* dMupx;
            vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:     ,:     ), N, 1), 0) .* dMupx;
            vecY1 = min(reshape(obj.U{ph,1+f}.y(:     ,1:Ny  ,:     ), N, 1), 0) .* dMupy;
            vecY2 = max(reshape(obj.U{ph,1+f}.y(:     ,2:Ny+1,:     ), N, 1), 0) .* dMupy;
            vecZ1 = min(reshape(obj.U{ph,1+f}.z(:     ,:     ,1:Nz  ), N, 1), 0) .* dMupz;
            vecZ2 = max(reshape(obj.U{ph,1+f}.z(:     ,:     ,2:Nz+1), N, 1), 0) .* dMupz;
            acc = Grid.Volume/dt .* (por .* obj.drhodp(Index.Start:Index.End) + rho .*dpor);
            
            DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J_PP = J_PP + spdiags(DiagVecs, DiagIndx, N, N);
            
            % 2. J_PT
            dMupx = obj.UpWind{ph,1+f}.x * ( obj.dMobdT(Index.Start:Index.End, ph) .* rho + obj.Mob(Index.Start:Index.End, ph) .* obj.drhodT(Index.Start:Index.End) );
            dMupy = obj.UpWind{ph,1+f}.y * ( obj.dMobdT(Index.Start:Index.End, ph) .* rho + obj.Mob(Index.Start:Index.End, ph) .* obj.drhodT(Index.Start:Index.End) );
            dMupz = obj.UpWind{ph,1+f}.z * ( obj.dMobdT(Index.Start:Index.End, ph) .* rho + obj.Mob(Index.Start:Index.End, ph) .* obj.drhodT(Index.Start:Index.End) );
            % Construct JPT block
            vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx  ,:     ,:     ), N, 1), 0) .* dMupx;
            vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:     ,:     ), N, 1), 0) .* dMupx;
            vecY1 = min(reshape(obj.U{ph,1+f}.y(:     ,1:Ny  ,:     ), N, 1), 0) .* dMupy;
            vecY2 = max(reshape(obj.U{ph,1+f}.y(:     ,2:Ny+1,:     ), N, 1), 0) .* dMupy;
            vecZ1 = min(reshape(obj.U{ph,1+f}.z(:     ,:     ,1:Nz  ), N, 1), 0) .* dMupz;
            vecZ2 = max(reshape(obj.U{ph,1+f}.z(:     ,:     ,2:Nz+1), N, 1), 0) .* dMupz;
            acc =  Grid.Volume/dt .* por .* obj.drhodT(Index.Start:Index.End) ;
            DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J_PT = spdiags(DiagVecs,DiagIndx,N,N);
            
            % Add Wells
            % for now, we will consider an only 2-phase system for adding the wells to the jacobian
            if f == 0 % only for reservoir
                [J_PP, J_PT, ~, ~] = obj.AddWellsToJacobian(J_PP, J_PT, Medium.State, Wells, Medium.K(:,1), ph);
            end
        end
        function [J_TP , J_TT] = BuildMediumHeatJacobian(obj, Medium, Wells, Grid, dt, Index, f, ph)
            % Create local variables
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N  = Grid.N;
            
            P   = Medium.State.Properties(['P_', num2str(ph)]).Value;
            T   = Medium.State.Properties(['T']).Value;
            rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
            h   = Medium.State.Properties(['h_', num2str(ph)]).Value;
            Rho_rock = Medium.Rho;

            Medium.ComputeDerPorosity(P);
            por = Medium.Por;
            dpor = Medium.DPor;
            
            % 1.J_TP
            % 1.a Pressure Block
            J_TP = obj.Th{ph, 1+f};
            
            % 1.b: compressibility part
            dMupx = obj.UpWind{ph,1+f}.x * ( obj.Mob(Index.Start:Index.End, ph) .* ( obj.drhodp(Index.Start:Index.End) .* h + obj.dhdp(Index.Start:Index.End) .* rho ) );
            dMupy = obj.UpWind{ph,1+f}.y * ( obj.Mob(Index.Start:Index.End, ph) .* ( obj.drhodp(Index.Start:Index.End) .* h + obj.dhdp(Index.Start:Index.End) .* rho ) );
            dMupz = obj.UpWind{ph,1+f}.z * ( obj.Mob(Index.Start:Index.End, ph) .* ( obj.drhodp(Index.Start:Index.End) .* h + obj.dhdp(Index.Start:Index.End) .* rho ) );
            
            vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx  ,:     ,:     ), N, 1), 0) .* dMupx;
            vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:     ,:     ), N, 1), 0) .* dMupx;
            vecY1 = min(reshape(obj.U{ph,1+f}.y(:     ,1:Ny  ,:     ), N, 1), 0) .* dMupy;
            vecY2 = max(reshape(obj.U{ph,1+f}.y(:     ,2:Ny+1,:     ), N, 1), 0) .* dMupy;
            vecZ1 = min(reshape(obj.U{ph,1+f}.z(:     ,:     ,1:Nz  ), N, 1), 0) .* dMupz;
            vecZ2 = max(reshape(obj.U{ph,1+f}.z(:     ,:     ,2:Nz+1), N, 1), 0) .* dMupz;
            acc = (Grid.Volume/dt) .* ( obj.Cp.*(por.*obj.drhodp(Index.Start:Index.End)+rho.*dpor) + Medium.Cpr.*(-dpor).*Rho_rock ) .* T;
            
            DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J_TP = J_TP + spdiags(DiagVecs, DiagIndx, N, N);
            
            % 2. J_TT
            J_TT = obj.Tk{1, 1+f};
            Mob  = obj.Mob(Index.Start:Index.End, ph);
            
            dMupx = obj.UpWind{ph,1+f}.x * (obj.dMobdT(Index.Start:Index.End) .* rho .* h + obj.drhodT(Index.Start:Index.End) .* Mob .* h  + obj.dhdT(Index.Start:Index.End) .* rho .* Mob);
            dMupy = obj.UpWind{ph,1+f}.y * (obj.dMobdT(Index.Start:Index.End) .* rho .* h + obj.drhodT(Index.Start:Index.End) .* Mob .* h  + obj.dhdT(Index.Start:Index.End) .* rho .* Mob);
            dMupz = obj.UpWind{ph,1+f}.z * (obj.dMobdT(Index.Start:Index.End) .* rho .* h + obj.drhodT(Index.Start:Index.End) .* Mob .* h  + obj.dhdT(Index.Start:Index.End) .* rho .* Mob);
            % Construct JTT block
            vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx  ,:     ,:     ), N, 1), 0) .* dMupx;
            vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:     ,:     ), N, 1), 0) .* dMupx;
            vecY1 = min(reshape(obj.U{ph,1+f}.y(:     ,1:Ny  ,:     ), N, 1), 0) .* dMupy;
            vecY2 = max(reshape(obj.U{ph,1+f}.y(:     ,2:Ny+1,:     ), N, 1), 0) .* dMupy;
            vecZ1 = min(reshape(obj.U{ph,1+f}.z(:     ,:     ,1:Nz  ), N, 1), 0) .* dMupz;
            vecZ2 = max(reshape(obj.U{ph,1+f}.z(:     ,:     ,2:Nz+1), N, 1), 0) .* dMupz;
            acc = (Grid.Volume/dt) .* ( obj.Cp .* por .*( obj.drhodT(Index.Start:Index.End) .* T + rho ) + Medium.Cpr.* (1-por) .* Rho_rock );
            
            DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J_TT = J_TT + spdiags(DiagVecs, DiagIndx, N, N);
            
            % Add Wells
            % for now, we will consider an only 2-phase system for adding the wells to the jacobian
            if f == 0 % only for reservoir
                [~, ~, J_TP, J_TT] = obj.AddWellsToJacobian(J_TP, J_TT, Medium.State, Wells, Medium.K(:,1), ph);
            end
        end
        function JacobianFull = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            %% Jacobian's assembly for single phase geothermal reservoir
            % | JPP  JPT | dP |    | RP |
            % | JTP  JTT | dT | = -| RT |
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            Reservoir = ProductionSystem.Reservoir;
            Fractures = ProductionSystem.FracturesNetwork.Fractures;
            Wells = ProductionSystem.Wells;       
            
            for ph=1:obj.NofPhases
                %% Jacobian of the reservoir 
                Index.Start = 1;
                Index.End = Nm;
                % Flow Jacobian blocks
                [J_PP, J_PT] = BuildMediumFlowJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0, ph); % 0 = reservoir
                % Heat Jacobian blocks
                [J_TP, J_TT] = BuildMediumHeatJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0, ph); % 0 = reservoir
                
                %% Jacobian of the fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Nf = DiscretizationModel.FracturesGrid.N;
                    Index.Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                    Index.End = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, Nf(f));
                    % Flow Jacobian Block
                    [Jf_PP, Jf_PT] = BuildMediumFlowJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f, ph);
                    J_PP  = blkdiag(J_PP, Jf_PP);
                    J_PT  = blkdiag(J_PT, Jf_PT);
                    % Heat Jacobian Block
                    [Jf_TP, Jf_TT] = BuildMediumHeatJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f, ph);
                    J_TP  = blkdiag(J_TP, Jf_TP);
                    J_TT  = blkdiag(J_TT, Jf_TT);
                end
                
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
                    T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                    T_Geo_Cond = DiscretizationModel.CrossConnections(c).T_Geo_Cond;
                    UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                    i = c + Nm;
                    j = DiscretizationModel.CrossConnections(c).Cells;
                    
                    %% J_PP Coupling
                    J_PP_1_conn = - T_Geo .*  UpWind(:, ph) .* obj.Mob(j, ph) .* ( rho(j,ph) + obj.drhodp(j, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_PP_2_conn = - T_Geo .* ~UpWind(:, ph) .* obj.Mob(i, ph) .* ( rho(i,ph) + obj.drhodp(i, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_PP_conn  = J_PP_1_conn + J_PP_2_conn;
                    
                    % frac - mat or frac1 - frac2
                    if sum(J_PP(i, j))~=0 || sum(J_PP(j, i))~=0
                        error('J_PP(i, j) or J_PP(j, i) is not zero!');
                    end
                    J_PP(i, j) = J_PP_conn;
                    J_PP(i, i) = J_PP(i, i) - sum(J_PP_conn);
                    % mat-frac or frac2 - frac1
                    J_PP(j, i) = J_PP_conn';
                    J_PP(sub2ind([Nt, Nt], j, j)) = J_PP(sub2ind([Nt, Nt], j, j)) - J_PP_conn;
                    
                    %% J_PT Coupling
                    J_PT_1_conn = - T_Geo .*  UpWind(:, ph) .* (P(j, ph) - P(i, ph)).* ( rho(j, ph) .* obj.dMobdT(j, ph) + obj.Mob(j, ph) .* obj.drhodT(j, ph) );
                    J_PT_2_conn = - T_Geo .* ~UpWind(:, ph) .* (P(j, ph) - P(i, ph)).* ( rho(i, ph) .* obj.dMobdT(i, ph) + obj.Mob(i, ph) .* obj.drhodT(i, ph) );
                    
                    % frac - mat or frac1 - frac2
                    if sum(J_PT(i, j))~=0 || sum(J_PT(j, i))~=0
                        error('J_PT(i, j) or J_PT(j, i) is not zero!');
                    end
                    J_PT(i, j) = J_PT_1_conn;
                    J_PT(i, i) = J_PT(i, i) + sum(J_PT_2_conn);
                    % mat-frac or frac2 - frac1
                    J_PT(j, i) = - J_PT_2_conn;
                    % diag of mat or frac2
                    J_PT(sub2ind([Nt, Nt], j, j)) = J_PT(sub2ind([Nt, Nt], j, j)) - J_PT_1_conn;
                    
                    %% J_TP Coupling
                    J_TP_1_conn = - T_Geo .*  UpWind(:, ph) .* obj.Mob(j, ph) .* ( (P(j, ph) - P(i, ph)) .* (rho(j, ph) .* obj.dhdp(j, ph) + h(j, ph) .* obj.drhodp(j, ph)) ...
                        + rho(j, ph) .* h(j, ph) );
                    J_TP_2_conn = - T_Geo .* ~UpWind(:, ph) .* obj.Mob(i, ph) .* ( (P(j, ph) - P(i, ph)) .* (rho(i, ph) .* obj.dhdp(i, ph) + h(i, ph) .* obj.drhodp(i, ph)) ...
                        + rho(i, ph) .* h(i, ph) );
                    J_TP_conn = J_TP_1_conn + J_TP_2_conn;
                    % frac - mat or frac1 - frac2
                    if sum(J_TP(i, j))~=0 || sum(J_TP(j, i))~=0
                        error('J_TP(i, j) or J_TP(j, i) is not zero!');
                    end
                    J_TP(i, j) = J_TP_conn;
                    J_TP(i, i) = J_TP(i, i) - sum(J_TP_conn);
                    % mat-frac or frac2 - frac1
                    J_TP(j, i) = J_TP_conn';
                    % diag of mat or frac2
                    J_TP(sub2ind([Nt, Nt], j, j)) = J_TP(sub2ind([Nt, Nt], j, j)) - J_TP_conn;
                    
                    %% J_TT Coupling
                    J_TT_1_conn = - T_Geo .*  UpWind(:, ph) .* ( (P(j, ph) - P(i, ph)) .* ( obj.dMobdT(j, ph) .* rho(j, ph) .* h(j, ph) + obj.Mob(j, ph) .* obj.drhodT(j, ph) .* h(j, ph) ...
                        + obj.Mob(j, ph) .* rho(j, ph) .* obj.dhdT(j, ph) ) ); - T_Geo_Cond .*  UpWind(:, ph);
                    J_TT_2_conn = - T_Geo .* ~UpWind(:, ph) .* ( (P(j, ph) - P(i, ph)) .* ( obj.dMobdT(i, ph) .* rho(i, ph) .* h(i, ph) + obj.Mob(i, ph) .* obj.drhodT(i, ph) .* h(i, ph) ...
                        + obj.Mob(i, ph) .* rho(i, ph) .* obj.dhdT(i, ph) ) ); - T_Geo_Cond .*  UpWind(:, ph);
                    % frac - mat or frac1 - frac2
                    if sum(J_TT(i, j))~=0 || sum(J_TT(j, i))~=0
                        error('J_TT(i, j) or J_TT(j, i) is not zero!');
                    end
                    J_TT(i, j) = J_TT_1_conn';
                    J_TT(i, i) = J_TT(i, i) + sum(J_TT_2_conn);
                    % mat-frac or frac2 - frac1
                    J_TT(j, i) = -J_TT_2_conn;
                    % diag of mat or frac2
                    J_TT(sub2ind([Nt, Nt], j, j)) = J_TT(sub2ind([Nt, Nt], j, j)) - J_TT_1_conn;
                end
            end
            
            % Build & Stack Jacobian
            JacobianFull = [J_PP, J_PT ; J_TP, J_TT];
        end
        function ConstrainedResidual(obj, FluidModel, ProductionSystem, DiscretizationModel, dt, State0)
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            Reservoir = ProductionSystem.Reservoir;
            Fractures = ProductionSystem.FracturesNetwork.Fractures;
            Wells = ProductionSystem.Wells;       
            ResidualFull = obj.BuildResidual(ProductionSystem, DiscretizationModel, dt, State0);
            
            %% Jacobian's assembly for single phase geothermal reservoir
            % | JPP  JPT | dP |    | RP |
            % | JTP  JTT | dT | = -| RT |
            for ph=1:obj.NofPhases
                %% Jacobian of the reservoir 
                Index.Start = 1;
                Index.End = Nm;
                % Flow Jacobian blocks
                [J_PP, ~] = BuildMediumFlowJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0, ph); % 0 = reservoir
                % Heat Jacobian blocks
                [~, J_TT] = BuildMediumHeatJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0, ph); % 0 = reservoir
                
                %% Jacobian of the fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Nf = DiscretizationModel.FracturesGrid.N;
                    Index.Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                    Index.End = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, Nf(f));
                    % Flow Jacobian Block
                    [Jf_PP, ~] = BuildMediumFlowJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f, ph);
                    J_PP  = blkdiag(J_PP, Jf_PP);
                    % Heat Jacobian Block
                    [~, Jf_TT] = BuildMediumHeatJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f, ph);
                    J_TT  = blkdiag(J_TT, Jf_TT);
                end
                
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
                    T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                    T_Geo_Cond = DiscretizationModel.CrossConnections(c).T_Geo_Cond;
                    UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                    i = c + Nm;
                    j = DiscretizationModel.CrossConnections(c).Cells;
                    
                    %% J_PP Coupling
                    J_PP_1_conn = - T_Geo .*  UpWind(:, ph) .* obj.Mob(j, ph) .* ( rho(j,ph) + obj.drhodp(j, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_PP_2_conn = - T_Geo .* ~UpWind(:, ph) .* obj.Mob(i, ph) .* ( rho(i,ph) + obj.drhodp(i, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_PP_conn  = J_PP_1_conn + J_PP_2_conn;
                    
                    % frac - mat or frac1 - frac2
                    if sum(J_PP(i, j))~=0 || sum(J_PP(j, i))~=0
                        error('J_PP(i, j) or J_PP(j, i) is not zero!');
                    end
                    J_PP(i, j) = J_PP_conn;
                    J_PP(i, i) = J_PP(i, i) - sum(J_PP_conn);
                    % mat-frac or frac2 - frac1
                    J_PP(j, i) = J_PP_conn';
                    J_PP(sub2ind([Nt, Nt], j, j)) = J_PP(sub2ind([Nt, Nt], j, j)) - J_PP_conn;
                    
                    %% J_TT Coupling
                    J_TT_1_conn = - T_Geo .*  UpWind(:, ph) .* ( (P(j, ph) - P(i, ph)) .* ( obj.dMobdT(j, ph) .* rho(j, ph) .* h(j, ph) + obj.Mob(j, ph) .* obj.drhodT(j, ph) .* h(j, ph) ...
                        + obj.Mob(j, ph) .* rho(j, ph) .* obj.dhdT(j, ph) ) ) - T_Geo_Cond .*  UpWind(:, ph);
                    J_TT_2_conn = - T_Geo .* ~UpWind(:, ph) .* ( (P(j, ph) - P(i, ph)) .* ( obj.dMobdT(i, ph) .* rho(i, ph) .* h(i, ph) + obj.Mob(i, ph) .* obj.drhodT(i, ph) .* h(i, ph) ...
                        + obj.Mob(i, ph) .* rho(i, ph) .* obj.dhdT(i, ph) ) ) - T_Geo_Cond .*  UpWind(:, ph);
                    % frac - mat or frac1 - frac2
                    if sum(J_TT(i, j))~=0 || sum(J_TT(j, i))~=0
                        error('J_TT(i, j) or J_TT(j, i) is not zero!');
                    end
                    J_TT(i, j) = J_TT_1_conn';
                    J_TT(i, i) = J_TT(i, i) + sum(J_TT_2_conn);
                    % mat-frac or frac2 - frac1
                    J_TT(j, i) = -J_TT_2_conn;
                    % diag of mat or frac2
                    J_TT(sub2ind([Nt, Nt], j, j)) = J_TT(sub2ind([Nt, Nt], j, j)) - J_TT_1_conn;
                end
            end
            
            %% Solve for delta
            J = blkdiag(J_PP,J_TT);
            delta = J\-ResidualFull;
            deltaP = delta(1:Nt);
            deltaT = delta(Nt+1:2*Nt);
            
            %% Update Reservoir State
            Pm = ProductionSystem.Reservoir.State.Properties('P_1');
            Pm.update(deltaP(1:Nm));
            Tm = ProductionSystem.Reservoir.State.Properties('T');
            Tm.update(deltaT(1:Nm));
            % Update Phase Densities
            FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
            % Update viscosity
            FluidModel.ComputePhaseViscosities(ProductionSystem.Reservoir.State);
            % Update enthalpy
            FluidModel.ComputePhaseEnthalpies(ProductionSystem.Reservoir.State);
            % Update wells
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
            
            %% Update fracture state
            % Update Pressure
            if ProductionSystem.FracturesNetwork.Active
                EP = Nm;
                Nf = DiscretizationModel.FracturesGrid.N;
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    IP = EP+1;
                    EP = IP + Nf(f) - 1;
                    Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('P_1');
                    Pf.update(deltaP(IP:EP));
                    Tf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('T');
                    Tf.update(deltaT(IP:EP));
                    
                    % Update remaining properties
                    % Update Phase Densities
                    FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    % Update viscosity
                    FluidModel.ComputePhaseViscosities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    % Update enthalpy
                    FluidModel.ComputePhaseEnthalpies(ProductionSystem.FracturesNetwork.Fractures(f).State);
                end
            end
        end
        function ConstrainedPressureResidual(obj, FluidModel, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute source terms
            [Qw, ~] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            Qf = zeros(DiscretizationModel.N, obj.NofPhases);
            if ProductionSystem.FracturesNetwork.Active
                [Qf, ~, ~] = obj.ComputeSourceTerms_frac_mat(ProductionSystem, DiscretizationModel);
            end
            % Initialise residual vector (Nph * N, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;
            if ProductionSystem.FracturesNetwork.Active
                Nf = DiscretizationModel.FracturesGrid.N;
            else
                Nf = 0;
            end
            Nt = DiscretizationModel.N;
            R_P = zeros( Nt , 1 );
            for ph=1:obj.NofPhases
                %% Transmissibilities of flow and heat
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                [obj.Tph{ph, 1}, obj.Gph{ph, 1}, ~, ~] = obj.TransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    ProductionSystem.Reservoir.State.Properties(['h_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
                % Fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End+1;
                    Index.End = Index.Start + Nf(f) - 1;
                    [obj.Tph{ph, 1+f}, obj.Gph{ph, 1+f},~, ~] = obj.TransmissibilityMatrix( ...
                        DiscretizationModel.FracturesGrid.Grids(f), ...
                        obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph), ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value, ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['h_',num2str(ph)]).Value, ...
                        obj.GravityModel.RhoInt{ph, 1+f});
                end
                
                %% Flow Residual
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                Residualm = BuildMediumFlowResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, Qw, Qf, 0, ph);
                R_P((ph-1)*Nt + Index.Start: (ph-1)*Nt + Index.End) = Residualm;
                % Fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End+1;
                    Index.End = Index.Start + Nf(f) - 1;
                    Residualf = BuildMediumFlowResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, Qw, Qf, f, ph);
                    R_P(Index.Start:Index.End) = Residualf;
                end
            end
            
            %% Jacobian's assembly for single phase geothermal reservoir
            % |JPP| |dP| = -|RP|
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            Reservoir = ProductionSystem.Reservoir;
            Fractures = ProductionSystem.FracturesNetwork.Fractures;
            Wells = ProductionSystem.Wells;
            
            for ph=1:obj.NofPhases
                %% Jacobian of the reservoir 
                Index.Start = 1;
                Index.End = Nm;
                % Flow Jacobian blocks
                [J_PP, ~] = BuildMediumFlowJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0, ph); % 0 = reservoir
                
                %% Jacobian of the fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Nf = DiscretizationModel.FracturesGrid.N;
                    Index.Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                    Index.End = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, Nf(f));
                    % Flow Jacobian Block
                    [Jf_PP, ~] = BuildMediumFlowJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f, ph);
                    J_PP  = blkdiag(J_PP, Jf_PP);
                end
                
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
                    
                    %% J_PP Coupling
                    J_PP_1_conn = - T_Geo .*  UpWind(:, ph) .* obj.Mob(j, ph) .* ( rho(j,ph) + obj.drhodp(j, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_PP_2_conn = - T_Geo .* ~UpWind(:, ph) .* obj.Mob(i, ph) .* ( rho(i,ph) + obj.drhodp(i, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_PP_conn  = J_PP_1_conn + J_PP_2_conn;
                    % frac - mat or frac1 - frac2
                    if sum(J_PP(i, j))~=0 || sum(J_PP(j, i))~=0
                        error('J_PP(i, j) or J_PP(j, i) is not zero!');
                    end
                    J_PP(i, j) = J_PP_conn;
                    J_PP(i, i) = J_PP(i, i) - sum(J_PP_conn);
                    % mat-frac or frac2 - frac1
                    J_PP(j, i) = J_PP_conn';
                    J_PP(sub2ind([Nt, Nt], j, j)) = J_PP(sub2ind([Nt, Nt], j, j)) - J_PP_conn;
                end
            end
            
            % Solve for deltaP
            if isprop(DiscretizationModel,'OperatorsHandler')
                deltaPc = ( DiscretizationModel.OperatorsHandler.ADMRest * J_PP * DiscretizationModel.OperatorsHandler.ADMProl{1} ) \ ...
                          (- DiscretizationModel.OperatorsHandler.ADMRest * R_P );
                deltaP = DiscretizationModel.OperatorsHandler.ADMProl{1} * deltaPc;
            else
                deltaP = J_PP\(-R_P);
            end
            
            %% Update Reservoir State
            Pm = ProductionSystem.Reservoir.State.Properties('P_1');
            Pm.update(deltaP(1:Nm));
            % Update Phase Densities
            FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
            % Update viscosity
            FluidModel.ComputePhaseViscosities(ProductionSystem.Reservoir.State);
            % Update enthalpy
            FluidModel.ComputePhaseEnthalpies(ProductionSystem.Reservoir.State);
            % Update wells
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
            
            %% Update fracture state
            % Update Pressure
            if ProductionSystem.FracturesNetwork.Active
                EP = Nm;
                Nf = DiscretizationModel.FracturesGrid.N;
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    IP = EP+1;
                    EP = IP + Nf(f) - 1;
                    Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_1']);
                    Pf.update(deltaP(IP:EP));
                    
                    % Update remaining properties
                    % Update Phase Densities
                    FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    % Update viscosity
                    FluidModel.ComputePhaseViscosities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    % Update enthalpy
                    FluidModel.ComputePhaseEnthalpies(ProductionSystem.FracturesNetwork.Fractures(f).State);
                end
            end
        end
        function ConstrainedTemperatureResidual(obj, FluidModel, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute source terms
            [~, Qhw] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            Qhf= zeros(DiscretizationModel.N, obj.NofPhases);              % Heat Convection flux betweem each two media
            RTf= zeros(DiscretizationModel.N, 1);                          % Heat Conduction flux betweem each two media
            if ProductionSystem.FracturesNetwork.Active
                [~, Qhf, RTf] = obj.ComputeSourceTerms_frac_mat(ProductionSystem, DiscretizationModel);
            end
            RTf= zeros(DiscretizationModel.N, 1);
            
            % Initialise residual vector (Nph * N, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;
            if ProductionSystem.FracturesNetwork.Active
                Nf = DiscretizationModel.FracturesGrid.N;
            else
                Nf = 0;
            end
            Nt = DiscretizationModel.N;
            R_T = zeros( Nt , 1 );
            for ph=1:obj.NofPhases
                %% Transmissibilities of flow and heat
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                [obj.Tph{ph, 1}, obj.Gph{ph, 1}, obj.Th{ph, 1}, obj.Tk{1, 1}] = obj.TransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    ProductionSystem.Reservoir.State.Properties(['h_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
                % Fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End+1;
                    Index.End = Index.Start + Nf(f) - 1;
                    [obj.Tph{ph, 1+f}, obj.Gph{ph, 1+f}, obj.Th{ph, 1+f}, obj.Tk{1, 1+f}] = obj.TransmissibilityMatrix( ...
                        DiscretizationModel.FracturesGrid.Grids(f), ...
                        obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph), ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value, ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['h_',num2str(ph)]).Value, ...
                        obj.GravityModel.RhoInt{ph, 1+f});
                end
                
                %% Heat Residual
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                Residualm = BuildMediumHeatResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, Qhw, Qhf, RTf, 0, ph);
                R_T(Index.Start: Index.End) = Residualm;
                % Fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End + 1;
                    Index.End = Index.Start + Nf(f) - 1;
                    Residualf = BuildMediumHeatResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, Qhw, Qhf, RTf, f, ph);
                    R_T(Index.Start:Index.End) = Residualf; 
                end
            end
            
            %% Jacobian's assembly for single phase geothermal reservoir
            % |JTT| |dT| = -|RT|
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            Reservoir = ProductionSystem.Reservoir;
            Fractures = ProductionSystem.FracturesNetwork.Fractures;
            Wells = ProductionSystem.Wells;
            
            for ph=1:obj.NofPhases
                %% Jacobian of the reservoir 
                Index.Start = 1;
                Index.End = Nm;
                % Heat Jacobian blocks
                [~, J_TT] = BuildMediumHeatJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0, ph); % 0 = reservoir
                
                %% Jacobian of the fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Nf = DiscretizationModel.FracturesGrid.N;
                    Index.Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                    Index.End = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, Nf(f));
                    % Heat Jacobian Block
                    [~, Jf_TT] = BuildMediumHeatJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f, ph);
                    J_TT  = blkdiag(J_TT, Jf_TT);
                end
                
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
                    T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                    T_Geo_Cond = DiscretizationModel.CrossConnections(c).T_Geo_Cond;
                    UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                    i = c + Nm;
                    j = DiscretizationModel.CrossConnections(c).Cells;
                    
                    %% J_TT Coupling
                    J_TT_1_conn = - T_Geo .*  UpWind(:, ph) .* ( (P(j, ph) - P(i, ph)) .* ( obj.dMobdT(j, ph) .* rho(j, ph) .* h(j, ph) + obj.Mob(j, ph) .* obj.drhodT(j, ph) .* h(j, ph) ...
                        + obj.Mob(j, ph) .* rho(j, ph) .* obj.dhdT(j, ph) ) ); - T_Geo_Cond .*  UpWind(:, ph);
                    J_TT_2_conn = - T_Geo .* ~UpWind(:, ph) .* ( (P(j, ph) - P(i, ph)) .* ( obj.dMobdT(i, ph) .* rho(i, ph) .* h(i, ph) + obj.Mob(i, ph) .* obj.drhodT(i, ph) .* h(i, ph) ...
                        + obj.Mob(i, ph) .* rho(i, ph) .* obj.dhdT(i, ph) ) ); - T_Geo_Cond .*  UpWind(:, ph);
                    % frac - mat or frac1 - frac2
                    if sum(J_TT(i, j))~=0 || sum(J_TT(j, i))~=0
                        error('J_TT(i, j) or J_TT(j, i) is not zero!');
                    end
                    J_TT(i, j) = J_TT_1_conn';
                    J_TT(i, i) = J_TT(i, i) + sum(J_TT_2_conn);
                    % mat-frac or frac2 - frac1
                    J_TT(j, i) = -J_TT_2_conn;
                    % diag of mat or frac2
                    J_TT(sub2ind([Nt, Nt], j, j)) = J_TT(sub2ind([Nt, Nt], j, j)) - J_TT_1_conn;
                end
            end
            
            % Solve for deltaT
            if isprop(DiscretizationModel,'OperatorsHandler')
                deltaTc = ( DiscretizationModel.OperatorsHandler.ADMRest * J_TT * DiscretizationModel.OperatorsHandler.ADMProl{2} ) \ ...
                          (- DiscretizationModel.OperatorsHandler.ADMRest * R_T );
                deltaT = DiscretizationModel.OperatorsHandler.ADMProl{2} * deltaTc;
            else
                deltaT = J_TT\(-R_T);
            end
            
            %% Update Reservoir State
            Tm = ProductionSystem.Reservoir.State.Properties('T');
            Tm.update(deltaT(1:Nm));
            % Update Phase Densities
            FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
            % Update viscosity
            FluidModel.ComputePhaseViscosities(ProductionSystem.Reservoir.State);
            % Update enthalpy
            FluidModel.ComputePhaseEnthalpies(ProductionSystem.Reservoir.State);
            % Update wells
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
            
            %% Update fracture state
            % Update Pressure
            if ProductionSystem.FracturesNetwork.Active
                EP = Nm;
                Nf = DiscretizationModel.FracturesGrid.N;
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    IP = EP+1;
                    EP = IP + Nf(f) - 1;
                    T = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('T');
                    T.update(deltaT(IP:EP));
                    
                    % Update remaining properties
                    % Update Phase Densities
                    FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    % Update viscosity
                    FluidModel.ComputePhaseViscosities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    % Update enthalpy
                    FluidModel.ComputePhaseEnthalpies(ProductionSystem.FracturesNetwork.Fractures(f).State);
                end
            end
        end
        function delta = UpdateState(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                return
            else
                delta = obj.UpdateTemperature(delta, ProductionSystem, FluidModel, DiscretizationModel);
                
                Nm = DiscretizationModel.ReservoirGrid.N;
                Nt = DiscretizationModel.N;
                deltaP = delta(1:Nt);
                deltaT = delta(Nt+1:2*Nt);
                
                %% Update reservoir state
                % 1. Update Pressure
                Pm = ProductionSystem.Reservoir.State.Properties(['P_1']);
                Pm.update(deltaP(1:Nm));
                % 2. Update Temperature
                Tm = ProductionSystem.Reservoir.State.Properties('T');
                Tm.update(deltaT(1:Nm));
                
                % Update remaining properties
                % 3. Update Phase Densities
                FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
                % 4. Update viscosity
                FluidModel.ComputePhaseViscosities(ProductionSystem.Reservoir.State);
                % 5. Update enthalpy
                FluidModel.ComputePhaseEnthalpies(ProductionSystem.Reservoir.State);
                % 6. Update wells
                ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
            
                %% Update fracture state
                % 1. Update Pressure
                if ProductionSystem.FracturesNetwork.Active
                    EP = Nm;
                    Nf = DiscretizationModel.FracturesGrid.N;
                    for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                        IP = EP+1;
                        EP = IP + Nf(f) - 1;
                        Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_1']);
                        Pf.update(deltaP(IP:EP));
                        % 2. Update Fluid Temperature
                        T = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('T');
                        T.update(deltaT(IP:EP));
                        
                        % Update remaining properties
                        % 3. Update Phase Densities
                        FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                        % 4. Update viscosity
                        FluidModel.ComputePhaseViscosities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                        % 5. Update enthalpy
                        FluidModel.ComputePhaseEnthalpies(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    end
                end
            end
        end
        function [Tph, Gph, Th, Tk] = TransmissibilityMatrix(obj, Grid, UpWind, Mob, rho, h, RhoInt)
            % Tph: Mass flow transmissibility
            % Gph: Gravity component of mass flow transmissibility
            % Th : Heat convection transmissibility (Tph*rho)
            % Tk : Heat conduction transmissibility
            
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            
            %% Tph matrix
            Tx = zeros(Nx+1, Ny, Nz);
            Ty = zeros(Nx, Ny+1, Nz);
            Tz = zeros(Nx, Ny, Nz+1);
            
            % Apply upwind operator
            Mupx = UpWind.x*(Mob .* rho);
            Mupy = UpWind.y*(Mob .* rho);
            Mupz = UpWind.z*(Mob .* rho);
            Mupx = reshape(Mupx, Nx, Ny, Nz);
            Mupy = reshape(Mupy, Nx, Ny, Nz);
            Mupz = reshape(Mupz, Nx, Ny, Nz);

            Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
            Tph = ReshapeTransmisibility(Grid, Nx, Ny, Nz, N, Tx, Ty, Tz); 
            
            %% Gph matrix
            Tx(2:Grid.Nx,:,:)= Tx(2:Grid.Nx,:,:) .* RhoInt.x(2:Grid.Nx,:,:);
            Ty(:,2:Grid.Ny,:)= Ty(:,2:Grid.Ny,:) .* RhoInt.y(:,2:Grid.Ny,:);
            Tz(:,:,2:Grid.Nz)= Tz(:,:,2:Grid.Nz) .* RhoInt.z(:,:,2:Grid.Nz);   
            Gph = ReshapeTransmisibility(Grid, Nx, Ny, Nz, N, Tx, Ty, Tz);
            
            %% Th Matrix
            Tx = zeros(Nx+1, Ny, Nz);
            Ty = zeros(Nx, Ny+1, Nz);
            Tz = zeros(Nx, Ny, Nz+1);

            % Apply upwind operator
            Mupx = UpWind.x*(Mob .* rho .* h);
            Mupy = UpWind.y*(Mob .* rho .* h);
            Mupz = UpWind.z*(Mob .* rho .* h);
            Mupx = reshape(Mupx, Nx, Ny, Nz);
            Mupy = reshape(Mupy, Nx, Ny, Nz);
            Mupz = reshape(Mupz, Nx, Ny, Nz);

            Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
            Th = ReshapeTransmisibility(Grid, Nx, Ny, Nz, N, Tx, Ty, Tz); 
            
            %% Tk Matrix
            THx = zeros(Nx+1, Ny, Nz);
            THy = zeros(Nx, Ny+1, Nz);
            THz = zeros(Nx, Ny, Nz+1);
            
            THx(2:Nx,:,:)= Grid.THx(2:Nx,:,:); 
            THy(:,2:Ny,:)= Grid.THy(:,2:Ny,:);
            THz(:,:,2:Nz)= Grid.THz(:,:,2:Nz);
            Tk = ReshapeTransmisibility(Grid, Nx, Ny, Nz, N, THx, THy, THz); % Transmisibility of rock conductivity

            % Construct matrix
            function Tph = ReshapeTransmisibility(Grid, Nx, Ny, Nz, N, Tx, Ty, Tz) % remove function - end
                x1 = reshape(Tx(1:Nx,:,:), N, 1);
                x2 = reshape(Tx(2:Nx+1,:,:), N, 1);
                y1 = reshape(Ty(:,1:Ny,:), N, 1);
                y2 = reshape(Ty(:,2:Ny+1,:), N, 1);
                z1 = reshape(Tz(:,:,1:Nz), N, 1);
                z2 = reshape(Tz(:,:,2:Nz+1), N, 1);
                DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
                DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
                Tph = spdiags(DiagVecs,DiagIndx,N,N);
            end
        end
        function [Qw, Qhw]= ComputeSourceTerms(obj, N, Wells)
            Qw = zeros(N, obj.NofPhases);
            Qhw = zeros(N, obj.NofPhases);  
            %Injectors
            for i=1:Wells.NofInj
                c = Wells.Inj(i).Cells;
                Qw(c, :) = Wells.Inj(i).QPhases(:,:);
                Qhw(c, :) = Wells.Inj(i).Qh(:,:);
            end
            %Producers
            for i=1:Wells.NofProd
                c = Wells.Prod(i).Cells;
                Qw(c, :) = Wells.Prod(i).QPhases(:,:);
                Qhw(c, :) = Wells.Prod(i).Qh(:,:);
            end
        end
        function [Qf, Qhf, RTf]= ComputeSourceTerms_frac_mat(obj, ProductionSystem, DiscretizationModel)
            Qf = zeros(DiscretizationModel.N, obj.NofPhases);              % Mass Flow flux between each two media 
            Qhf = zeros(DiscretizationModel.N, obj.NofPhases);             % Heat Convection flux betweem each two media
            RTf = zeros(DiscretizationModel.N, 1);                         % Heat Conduction flux betweem each two media
            Nm = DiscretizationModel.ReservoirGrid.N;

            % Global variables
            P = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'P_'); % useful for cross connections assembly
            rho = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'rho_'); % useful for cross connections assembly
            h = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'h_'); % useful for cross connections assembly
            T = ProductionSystem.CreateGlobalSinglePhaseVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], 'T'); % useful for cross connections assembly
            for c=1:length(DiscretizationModel.CrossConnections)
                j = DiscretizationModel.CrossConnections(c).Cells;
                i = c + Nm;
                T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                T_Geo_Cond = DiscretizationModel.CrossConnections(c).T_Geo_Cond;
                UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                for ph=1:obj.NofPhases
                    Qf(i, ph) = Qf(i, ph) + sum(T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* (P(j, ph) - P(i, ph)) + ...
                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* (P(j, ph) - P(i, ph))) );
                    Qf(j, ph) = Qf(j, ph) + T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* (P(i, ph) - P(j, ph)) + ...
                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* (P(i, ph) - P(j, ph)) );
                    
                    Qhf(i, ph) = Qhf(i, ph) + sum(T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* h(j,ph) .* (P(j, ph) - P(i, ph)) + ...
                        ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* h(i,ph) .* (P(j, ph) - P(i, ph))) );
                    Qhf(j, ph) = Qhf(j, ph) + T_Geo .* ...
                        ( UpWind(:, ph) .* obj.Mob(j, ph) .* rho(j, ph) .* h(j,ph) .* (P(i, ph) - P(j, ph)) + ...
                         ~UpWind(:, ph) .* obj.Mob(i, ph) .* rho(i, ph) .* h(i,ph) .* (P(i, ph) - P(j, ph)) );
                end
                RTf(i) = RTf(i) + sum( T_Geo_Cond .* (T(j) - T(i)) );
                RTf(j) = RTf(j) + T_Geo_Cond .* (T(i) - T(j));
            end
        end
        function [AP1, AT1, AP2, AT2] = AddWellsToJacobian(obj, AP, AT, State, Wells, K, ph)
            % Define Local handles
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            AP1 = AP; % add wells to pressure-pressure jacobian (mass balance)
            AT1 = AT; % add wells to pressure-heat jacobian (mass balance)
            AP2 = AP; % add wells to heat-pressure jacobian (energy balance)
            AT2 = AT; % add wells to heat-heat jacobian (energy balance)
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                [dQdp, ~] = Inj(i).dQdPdT(K, obj.NofPhases);
                [dQhdp, ~] = Inj(i).dQhdPdT(K, obj.NofPhases);
                for j=1:length(a)
                    % add derivative of inj well to Jpp
                    AP1(a(j),a(j)) = AP1(a(j),a(j)) - dQdp(j, ph);
                    % add derivative of inj well to JTp
                    AP2(a(j),a(j)) = AP2(a(j),a(j)) - dQhdp(j, ph);
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                [dQdp , dQdT ] = Prod(i).dQdPdT(State, K, obj.Mob, obj.dMobdT, obj.drhodp, obj.drhodT, obj.NofPhases);
                [dQhdp, dQhdT] = Prod(i).dQhdPdT(State, K, obj.Mob, obj.dMobdT, obj.drhodp, obj.drhodT, obj.dhdp, obj.dhdT, obj.NofPhases);
                for j=1:length(b)
                    % add derivative of prod well to Jpp & JpT
                    AP1(b(j),b(j)) = AP1(b(j),b(j)) - dQdp(j, ph);                    
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