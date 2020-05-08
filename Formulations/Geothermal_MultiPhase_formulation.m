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
        MatrixAssembler
        
        drhoTdp
        drho_times_Sdp
        drho_times_hdp
        dMobdp
        dSdp
        drhoHSdp
        dTdp

        drhodh
        drhoTdh
        drho_times_Sdh
        drho_times_hdh
        dMobdh
        dSdh
        drhoHSdh
        d2rhodh2
        d2Mobdh2
        dTdh
        
%         Kf % fluid thermal conductivity
        Tk % transmisibility of rock conductivity
        Th % transmisibility of rho .* h
        
        Thph
        Ghph
    end
    methods
        function obj = Geothermal_MultiPhase_formulation(NofPhases)
            obj@formulation();
            obj.Tph = cell(NofPhases,1);
            obj.Gph = cell(NofPhases,1);
            obj.MatrixAssembler = matrix_assembler_geothermal();
            % This is the matrix assembler for geothermal (multiphase)
        end
        function x = GetPrimaryUnknowns(obj, ProductionSystem, DiscretizationModel)
             Nt = DiscretizationModel.N; % total grid
             Nm = DiscretizationModel.ReservoirGrid.N; % matrix grid
             x = zeros(2 * Nt, 1); % 2 unknowns: P, H
             Index.Start = 1;
             Index.End = Nm;
                     
             % Get matrix unknowns
             x(Index.Start     :Index.End       ) = ProductionSystem.Reservoir.State.Properties('P_1').Value;
             x(Nt+Index.Start  :Nt + Index.End  ) = ProductionSystem.Reservoir.State.Properties('hTfluid').Value;
             
             % Get fracture unknowns
             for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                 Index.Start = Index.End+1;
                 Index.End = Index.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                 x(Index.Start     :Index.End       ) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('P_1').Value;
                 x(Nt+Index.Start  :Nt + Index.End  ) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('hTfluid').Value;
             end
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            % 1. Computing properties
            obj.ComputeProperties(ProductionSystem, FluidModel);
            
            % 2. Reservoir Properties and Derivatives
            obj.ComputeDerivatives(ProductionSystem, FluidModel);
        end
        function ComputeProperties(obj, ProductionSystem, FluidModel)
            %% 1. Geothermal Properties

            %% 2. Reservoir Properties and Derivatives
            % Reservoir
            FluidModel.GetTableIndex(ProductionSystem.Reservoir.State); % Get table indices for reservoir
            FluidModel.ComputeTableGrid();
            FluidModel.GetPHValues(ProductionSystem.Reservoir.State);
            
            FluidModel.GetPhaseDensities(ProductionSystem.Reservoir.State); 
            FluidModel.GetPhaseSaturations(ProductionSystem.Reservoir.State); 
            FluidModel.GetPhaseViscosities(ProductionSystem.Reservoir.State); 
            FluidModel.GetPhaseThermalConductivities(ProductionSystem.Reservoir.State);
            FluidModel.GetPhaseEnthalpies(ProductionSystem.Reservoir.State);
            
            FluidModel.GetTotalDensity(ProductionSystem.Reservoir.State);
            FluidModel.GetTemperature(ProductionSystem.Reservoir.State);
            
            % Compute enthalpy of rock
            FluidModel.ComputeRockEnthalpy(ProductionSystem.Reservoir);
            % Compute Effective Thermal Conductivity coefficient
            FluidModel.ComputeThermalConductivity(ProductionSystem.Reservoir);
            
            % for phase mobilities; including rel.perm model (linear). Computed based on single saturation value as S = Sa+Sb
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State);
                        
            % Update Pc [THIS GUY UPDATES P_1; Pfffff....]
            FluidModel.ComputePc(ProductionSystem.Reservoir.State);
            
            %% 3. Fractures Properties
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
%                 obj.drhodp = [obj.drhodp; FluidModel.ComputeDrhoDp(ProductionSystem.FracturesNetwork.Fractures(f).State)];
            end
        end
        function ComputeDerivatives(obj, ProductionSystem, FluidModel)
            % Pressure and Enthalpy indices for tables are updated when
            % updating the properties; if properties are not updated, there
            % is no point updating the derivatives
            
            % Derivatives with respect to Pressure
            obj.drhodp = FluidModel.ComputeDrhoDp(); 
            obj.drhoTdp = FluidModel.ComputeDrhoTDp(); 
            obj.drho_times_Sdp = FluidModel.ComputeDrho_times_SDp(); 
            obj.drho_times_hdp = FluidModel.ComputeDrho_times_hDp();   
            obj.dMobdp = FluidModel.ComputeDMobDp(ProductionSystem.Reservoir.State);
            obj.dSdp = FluidModel.ComputeDSDp();
            obj.drhoHSdp = FluidModel.ComputeDrhoHSDp();
            obj.dTdp = FluidModel.ComputeDTDp();
            
            % Derivatives with respect to Enthalpy
            obj.drhodh = FluidModel.ComputeDrhoDh(); 
            obj.drhoTdh = FluidModel.ComputeDrhoTDh(); 
            obj.drho_times_Sdh = FluidModel.ComputeDrho_times_SDh();
            obj.drho_times_hdh = FluidModel.ComputeDrho_times_hDh();   
            obj.dMobdh = FluidModel.ComputeDMobDh(ProductionSystem.Reservoir.State);
            obj.dSdh = FluidModel.ComputeDSDh();
            obj.drhoHSdh = FluidModel.ComputeDrhoHSDh();
            obj.dTdh = FluidModel.ComputeDTDh();
            
            
            %% 3. Fractures Derivatives
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
%                 obj.drhodp = [obj.drhodp; FluidModel.ComputeDrhoDp(ProductionSystem.FracturesNetwork.Fractures(f).State)];
            end
        end
        function Residual_MB = BuildMediumFlowResidual(obj, Medium, Grid, dt, State0, Index, qw, qf, f)
            depth = Grid.Depth;
            
            % Initialization with zeros
            Residual_MB = zeros(Grid.N,1);
            
            for ph=1:obj.NofPhases
                % Create local variables
                % here, you are using state0 for the old pressure so makes
                % sense
                P_old = State0.Properties(['P_', num2str(ph)]).Value(Index.Start:Index.End);                
                rho_old = State0.Properties(['rho_', num2str(ph)]).Value(Index.Start:Index.End);
                S_old = State0.Properties(['S_', num2str(ph)]).Value(Index.Start:Index.End);

                % here, you get the actual(updated) properties at the
                % current time step
                P_new = Medium.State.Properties(['P_', num2str(ph)]).Value;
                rho_new = Medium.State.Properties(['rho_', num2str(ph)]).Value;
                S_new = Medium.State.Properties(['S_', num2str(ph)]).Value;

                Medium.ComputePorosity(P_old);
                pv_old = Medium.Por*Grid.Volume;
                Medium.ComputePorosity(P_new);
                pv_new = Medium.Por*Grid.Volume;

                % Accumulation Term
                Accumulation = (pv_new.*rho_new.*S_new - pv_old.*rho_old.*S_old)/dt;
                % this could also be written using rhoT
                
                
                % RESIDUAL ( = LHS - RHS )               
                Residual_MB  = Residual_MB + Accumulation ...
                    + obj.Tph{ph, 1+f} * P_new...
                    - obj.Gph{ph, 1+f} * depth...
                    - qw(Index.Start:Index.End, ph)...
                    - qf(Index.Start:Index.End, ph);
            end
        end
        function Residual_EB = BuildMediumHeatResidual(obj, Medium, Grid, dt, State0, Index, qhw, qhf, RTf, f)
            T_old = State0.Properties('T').Value(Index.Start:Index.End); % T at previous time step
            T_new = Medium.State.Properties('T').Value;
            Rho_rock = Medium.Rho;
            depth = Grid.Depth;
            
            P_old = State0.Properties(['P_', num2str(1)]).Value(Index.Start:Index.End); 
            P_new = Medium.State.Properties(['P_', num2str(1)]).Value;
            H_new = Medium.State.Properties('hTfluid').Value(Index.Start:Index.End);
            % Above; a small trick is applied to overcome the issue with
            % phase pressure implementation, i.e. mv_(old,new), in the rock accumulation.
            % --> This will be resolved when residuals for each phase are
            % added together. Right now, it means capillary pressure is not
            % an option.

            % Pore Volume & Rock Volume
            Medium.ComputePorosity(P_old);
            mv_old = (1 - Medium.Por) * Grid.Volume; % Old rock volume
            Medium.ComputePorosity(P_new);           % Updating porosity
            mv_new = (1 - Medium.Por) * Grid.Volume; % New rock volume
            
            % Accumulation rock
            h_rock_old = mv_old .* Rho_rock .* Medium.Cpr .* T_old; % this can be written with hRock property
            h_rock_new = mv_new .* Rho_rock .* Medium.Cpr .* T_new;
            Accumulation_rock = (h_rock_new - h_rock_old)/dt;

            % Initialization with rock accumulation
            Residual_EB = Accumulation_rock;
%             Residual_EB = Accumulation_rock + obj.Tk{1, 1+f} * T_new - RTf(Index.Start:Index.End, 1);
%             Residual_EB = Accumulation_rock + obj.Tk{1, 1+f} * obj.dTdp(:,1) .* P_new + obj.Tk{1, 1+f} * obj.dTdh .* H_new;

            % Above; we add terms to the residual that are independent of
            % phase, i.e. rock-enthalpy and thermal conductivity (this value is the total conductivity, i.e. rock+water+steam)
            % --> : I LOOP OVER THE PHASES IN THE FUNCTION THAT COMPUTES THERM.COND TENSOR D !! i.e. effective thermal conductivity
            % (1-phi)*D_rock + phi*S_water*D_water + phi*S_steam*D_steam
            
            
            for ph=1:obj.NofPhases
                % Create local variables
                rho_old = State0.Properties(['rho_', num2str(ph)]).Value(Index.Start:Index.End);
                P_old = State0.Properties(['P_', num2str(ph)]).Value(Index.Start:Index.End); %Medium.State.Properties(['P_', num2str(ph)]).Value;
                S_old = State0.Properties(['S_', num2str(ph)]).Value(Index.Start:Index.End);
                h_old = State0.Properties(['h_', num2str(ph)]).Value(Index.Start:Index.End);

                rho_new = Medium.State.Properties(['rho_', num2str(ph)]).Value;
                P_new = Medium.State.Properties(['P_', num2str(ph)]).Value;
                S_new = Medium.State.Properties(['S_', num2str(ph)]).Value;
                h_new = Medium.State.Properties(['h_', num2str(ph)]).Value;

                % This part incorporates phase pressures, so theoretically
                % this does apply to capillary pressure
                Medium.ComputePorosity(P_old);
                pv_old = Medium.Por*Grid.Volume;         % Old pore Volume
                Medium.ComputePorosity(P_new);           % Updating porosity
                pv_new = Medium.Por*Grid.Volume;         % New pore volume
                
                % Accumulation fluid
                h_fluid_old = pv_old .* rho_old .* h_old .* S_old;
                h_fluid_new = pv_new .* rho_new .* h_new .* S_new;
                Accumulation_fluid = (h_fluid_new - h_fluid_old)/dt;
                % WHY DONT WE USE HTFLUID HERE ??? --> SHOULD NOT MAKE A
                % DIFFERENCE, BUT WHO KNOWS, MAYBE IT DOES...?

                % RESIDUAL
                Residual_EB  = Residual_EB + Accumulation_fluid ...
                    + obj.Thph{ph, 1+f} * P_new ...
                    - obj.Ghph{ph, 1+f} * depth ...
                    - qhw(Index.Start:Index.End, ph)...
                    - qhf(Index.Start:Index.End, ph);
            end
        end
        function ResidualFull = BuildResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute source terms
            [Qw, Qhw] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            Qf = zeros(DiscretizationModel.N, obj.NofPhases);              % Mass Flow flux between each two media
            Qhf= zeros(DiscretizationModel.N, obj.NofPhases);              % Heat Convection flux betweem each two media
            RTf= zeros(DiscretizationModel.N, 1);                          % Heat Conduction flux betweem each two media
            
            % Initialise residual vector (2 * N, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            ResidualFull = zeros( 2*Nt , 1 ); % We have only 2 equations now
            
            % Compute conductive tranmissibility
            obj.Tk{1,1} = obj.MatrixAssembler.ConductiveHeatTransmissibilityMatrix( DiscretizationModel.ReservoirGrid , ProductionSystem.Reservoir.State );
            
            for ph=1:obj.NofPhases
                %% Transmissibilities of flow and heat
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                % Compute flow transmissibility and gravity
                [obj.Tph{ph, 1}, obj.Gph{ph, 1}] = obj.MatrixAssembler.TransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
                % Compute (convective) heat transmissibility and gravity
                [obj.Thph{ph, 1}, obj.Ghph{ph, 1}] = obj.MatrixAssembler.ConvectiveHeatTransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    ProductionSystem.Reservoir.State.Properties(['h_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
            end
            %% Flow Residual
            % Reservoir
            Index.Start = 1;
            Index.End = Nm;
            Residualm = BuildMediumFlowResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, Qw, Qf, 0);
            ResidualFull(Index.Start: Index.End) = Residualm;
            
            %% Heat Residual
            % Reservoir
            Index.Start = Index.End + 1;
            Index.End = Index.Start + Nm - 1;
            Index_r.Start = 1;
            Index_r.End = Nm;
            Residualm = BuildMediumHeatResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index_r, Qhw, Qhf, RTf, 0);
            ResidualFull(Index.Start: Index.End) = Residualm;
        end
        function [J_MB_P , J_MB_H] = BuildMediumFlowJacobian(obj, Medium, Wells, Grid, dt, Index, f)
            % Create local variables
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            
            J_MB_P = sparse(N,N);
            J_MB_H = sparse(N,N);
            
            for ph=1:obj.NofPhases
                
                P = Medium.State.Properties(['P_', num2str(ph)]).Value;
                rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
                S = Medium.State.Properties(['S_',num2str(ph)]).Value;
                Medium.ComputeDerPorosity(P);
                phi = Medium.Por;
                dphidp = Medium.DPor;
                            
                %% 1.J_MB_P Block
                % 1.a Transmissibility part
                J_MB_P = J_MB_P + obj.Tph{ph,1+f};
                
                % 1.b: Accumulation part
                acc = (Grid.Volume/dt) .* (phi .* obj.drho_times_Sdp(Index.Start:Index.End,ph) + rho .* S .*dphidp);
                
                % 1.c: Additional derivative term (P^nu) in convective flux
                dMupx = obj.UpWind{ph,1+f}.x *( obj.Mob(Index.Start:Index.End, ph)    .* obj.drhodp(Index.Start:Index.End, ph) + ...
                                                obj.dMobdp(Index.Start:Index.End, ph) .* rho(Index.Start:Index.End)            );
                dMupy = obj.UpWind{ph,1+f}.y *( obj.Mob(Index.Start:Index.End, ph)    .* obj.drhodp(Index.Start:Index.End, ph) + ...
                                                obj.dMobdp(Index.Start:Index.End, ph) .* rho(Index.Start:Index.End)            );
                dMupz = obj.UpWind{ph,1+f}.z *( obj.Mob(Index.Start:Index.End, ph)    .* obj.drhodp(Index.Start:Index.End, ph) + ...
                                                obj.dMobdp(Index.Start:Index.End, ph) .* rho(Index.Start:Index.End)            );
                
                % Because of multiplication with velocity obj.U, we are multiplying it with grad(P^nu).
                vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                vecY1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                vecY2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                vecZ1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                vecZ2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
                % zeros(N,1); %

                % construction of J_MB_P
                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                J_MB_P = J_MB_P + spdiags(DiagVecs, DiagIndx, N, N);
                
                
                %% 2. J_MB_H Block
                % 2.a: Transmissibility part (zero because d(P^nu)/dH = 0)
                J_MB_H = J_MB_H + 0;
                
                % 2.b: Accumulation part
                acc = (Grid.Volume/dt) .* phi .* obj.drho_times_Sdh(Index.Start:Index.End,ph);
                
                % 2.c: Additional derivative term (P^nu) in convective flux; we get d(rho*lambda)/dH * P^nu
                dMupx = obj.UpWind{ph,1+f}.x * ( obj.Mob(Index.Start:Index.End, ph)    .* obj.drhodh(Index.Start:Index.End, ph) + ...
                                                 obj.dMobdh(Index.Start:Index.End, ph) .* rho(Index.Start:Index.End)            );
                dMupy = obj.UpWind{ph,1+f}.y * ( obj.Mob(Index.Start:Index.End, ph)    .* obj.drhodh(Index.Start:Index.End, ph) + ...
                                                 obj.dMobdh(Index.Start:Index.End, ph) .* rho(Index.Start:Index.End)            );
                dMupz = obj.UpWind{ph,1+f}.z * ( obj.Mob(Index.Start:Index.End, ph)    .* obj.drhodh(Index.Start:Index.End, ph) + ...
                                                 obj.dMobdh(Index.Start:Index.End, ph) .* rho(Index.Start:Index.End)            );
%                 dMupx = zeros(N,1);
%                 dMupy = zeros(N,1);
%                 dMupz = zeros(N,1);
                
                % Because of multiplication with velocity obj.U, we are multiplying it with grad(P^nu).
                vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                vecY1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                vecY2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                vecZ1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                vecZ2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
                
                % construction of J_MB_H
                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                J_MB_H = J_MB_H + spdiags(DiagVecs, DiagIndx, N, N);
            end
                            
            % Add Wells
            % for now, we will consider an only 2-phase system for adding the wells to the jacobian
            if f == 0 % only for reservoir
                [J_MB_P, J_MB_H] = obj.AddWellsToMassBalanceJacobian(J_MB_P, J_MB_H, Medium.State, Wells, Medium.K(:,1));
            end
        end
        function [J_EB_P , J_EB_H] = BuildMediumHeatJacobian(obj, Medium, Wells, Grid, dt, Index, f)
            % Create local variables
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N  = Grid.N;
            
            J_EB_P = sparse(N,N);
            J_EB_H = sparse(N,N);
            
            Rho_rock = Medium.Rho;
            phi = Medium.Por;
            Medium.ComputeDerPorosity(Medium.State.Properties('P_2').Value);
            dphidp = Medium.DPor;
                
            %% 1. J_EB_P Block
            % 1.a: Accumulation part (only for rock)
            acc_rock = (Grid.Volume/dt) .* (-1) .* dphidp .* Rho_rock .* Medium.State.Properties('hRock').Value; 
            J_EB_P = J_EB_P + spdiags(acc_rock, 0, N, N);

            % Heat conduction flux
%             J_EB_P = J_EB_P + obj.Tk{1,1+f} * obj.dTdp;

            % 1.b: Heat convection flux
            for ph=1:obj.NofPhases
                % Original transmissibility part
                J_EB_P = J_EB_P + obj.Thph{ph, 1+f};
                
                rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
                h = Medium.State.Properties(['h_',num2str(ph)]).Value;
                S = Medium.State.Properties(['S_',num2str(ph)]).Value;
                
                % 1.c: Accumulation part (only for fluid)
                acc_fluid = (Grid.Volume/dt) .* ( dphidp .* rho .* h .* S + phi .* obj.drhoHSdp(:,ph) );
                
                % 1.d: Additional derivative term (P^nu) in convective flux
                dMupx = obj.UpWind{ph,1+f}.x*( obj.Mob(Index.Start:Index.End, ph) .* obj.drho_times_hdp(Index.Start:Index.End, ph) + ...
                    obj.dMobdp(Index.Start:Index.End, ph) .* rho(Index.Start:Index.End) .* h(Index.Start:Index.End) );
                dMupy = obj.UpWind{ph,1+f}.y*( obj.Mob(Index.Start:Index.End, ph) .* obj.drho_times_hdp(Index.Start:Index.End, ph) + ...
                    obj.dMobdp(Index.Start:Index.End, ph) .* rho(Index.Start:Index.End) .* h(Index.Start:Index.End) );
                dMupz = obj.UpWind{ph,1+f}.z*( obj.Mob(Index.Start:Index.End, ph) .* obj.drho_times_hdp(Index.Start:Index.End, ph) + ...
                    obj.dMobdp(Index.Start:Index.End, ph) .* rho(Index.Start:Index.End) .* h(Index.Start:Index.End) );
                
                % Because of multiplication with velocity obj.U, we are multiplying it with grad(P^nu).
                vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                vecY1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                vecY2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                vecZ1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                vecZ2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
                % zeros(N,1); %

                % construction of J_MB_P
                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc_fluid, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                J_EB_P = J_EB_P + spdiags(DiagVecs, DiagIndx, N, N);
            end

                        
            %% 2. J_EB_H Block
            % 2.a Accumulation fluid part 1
            rhoT = Medium.State.Properties('rhoT').Value;
            acc_fluid = phi .* rhoT .* (Grid.Volume/dt);
            J_EB_H = J_EB_H + spdiags(acc_fluid, 0, N, N);
            
            % Heat conduction flux
%             J_EB_P = J_EB_P + obj.Tk{1,1+f} * obj.dTdh;
           
            % 2.b Heat convection and conduction fluxes
            for ph=1:obj.NofPhases
                % Original transmissibility part                
                J_EB_H = J_EB_H + 0;
                
                rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
                h = Medium.State.Properties(['h_',num2str(ph)]).Value;
                S = Medium.State.Properties(['S_',num2str(ph)]).Value;
                hT = Medium.State.Properties('hTfluid').Value;
                
                % 2.c: Accumulation fluid part 2 
%                 acc_fluid = (Grid.Volume/dt) .* ( phi .* obj.drhoHSdh(:,ph) );
                % if you multiply tables with phase enthalpy, the phase enthalpy becomes a function of 
                % enthalpy as well. That is what we did... However, it is
                % not helping the convergence, at all !
                
                % it should be:
                acc_fluid = zeros(N,1); 
                
                % Secondary option: THIS ONE ALSO WORKS !!
%                 acc_fluid = (Grid.Volume/dt) .* ( phi .* obj.drho_times_Sdh(:,ph) .* hT );
                  % --> So that means, times  rho .* S .* h ./ rhoT; 
                  % where you should sum over the derivative, but divide
                  % by rhoT only once...

                % 2.d: Additional derivative term (P^nu) in convective flux
                dMupx = obj.UpWind{ph,1+f}.x*( obj.Mob(Index.Start:Index.End, ph) .* obj.drhodh(Index.Start:Index.End, ph) + ...
                    obj.dMobdh(Index.Start:Index.End, ph) .* rho(Index.Start:Index.End) ) .* h(Index.Start:Index.End);
                dMupy = obj.UpWind{ph,1+f}.y*( obj.Mob(Index.Start:Index.End, ph) .* obj.drhodh(Index.Start:Index.End, ph) + ...
                    obj.dMobdh(Index.Start:Index.End, ph) .* rho(Index.Start:Index.End) ) .* h(Index.Start:Index.End);
                dMupz = obj.UpWind{ph,1+f}.z*( obj.Mob(Index.Start:Index.End, ph) .* obj.drhodh(Index.Start:Index.End, ph) + ...
                    obj.dMobdh(Index.Start:Index.End, ph) .* rho(Index.Start:Index.End) ) .* h(Index.Start:Index.End);
%                 dMupx = zeros(N,1);
%                 dMupy = zeros(N,1);
%                 dMupz = zeros(N,1);

                % Because of multiplication with velocity obj.U, we are multiplying it with grad(P^nu).
                vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                vecY1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                vecY2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                vecZ1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                vecZ2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
                % zeros(N,1); %

                % construction of J_MB_P
                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc_fluid, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                J_EB_H = J_EB_H + spdiags(DiagVecs, DiagIndx, N, N);
            end

            %% Add Wells
            % for now, we will consider an only 2-phase system for adding the wells to the jacobian
            if f == 0 % only for reservoir
                [J_EB_P, J_EB_H] = obj.AddWellsToEnergyBalanceJacobian(J_EB_P, J_EB_H, Medium.State, Wells, Medium.K(:,1));
            end
        end
        function JacobianFull = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            %% Jacobian's assembly for MultiPhase geothermal reservoir
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
            
            %% Jacobian of the reservoir
            Index.Start = 1;
            Index.End = Nm;
            % Flow Jacobian blocks
            [J_MB_P, J_MB_H] = BuildMediumFlowJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0); % 0 = reservoir
            % Heat Jacobian blocks
            [J_EB_P, J_EB_H] = BuildMediumHeatJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0); % 0 = reservoir
            
            % Build & Stack Jacobian
            JacobianFull = [J_MB_P, J_MB_H ; J_EB_P, J_EB_H];
        end
        function ConstrainedPressureResidual(obj, FluidModel, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute source terms
            [Qw, ~] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            Qf = zeros(DiscretizationModel.N, obj.NofPhases);
            % Initialise residual vector (Nph * N, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            for ph=1:obj.NofPhases
                %% Transmissibilities of flow and heat
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                % Compute flow transmissibility and gravity
                [obj.Tph{ph, 1}, obj.Gph{ph, 1}] = obj.MatrixAssembler.TransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
            end
                
            %% Mass Balance Residual
            Index.Start = 1;
            Index.End = Nm;
            R_MB = BuildMediumFlowResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, Qw, Qf, 0);
            
            %% Mass Balance Jacobian (J_MB_P)
            Nm = DiscretizationModel.ReservoirGrid.N;
            Reservoir = ProductionSystem.Reservoir;
            Wells = ProductionSystem.Wells;
            Index.Start = 1;
            Index.End = Nm;
            % Flow Jacobian blocks
            [J_MB_P, ~] = BuildMediumFlowJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0); % 0 = reservoir

            % Solve for deltaP
            deltaP = J_MB_P\(-R_MB);
                        
            %% Update Reservoir State
            Pm = ProductionSystem.Reservoir.State.Properties('P_2');
            Pm.update(deltaP(1:Nm));
            % Update ALL properties
            obj.ComputeProperties(ProductionSystem, FluidModel)            
            % 7. Update wells
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
            
            fprintf('Norm CPR = %1.5e \n', norm(R_MB) );
            fprintf('Norm Delta P = %1.5e \n', norm(deltaP./max(Pm.Value)) );

        end
        function ConstrainedEnthalpyResidual(obj, FluidModel, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute source terms
            [~,Qhw] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            Qhf = zeros(DiscretizationModel.N, obj.NofPhases);
            RTf = zeros(DiscretizationModel.N, obj.NofPhases);
            % Initialise residual vector (Nph * N, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;
            Nt = DiscretizationModel.N;
            % Compute conductive tranmissibility
            obj.Tk{1,1} = obj.MatrixAssembler.ConductiveHeatTransmissibilityMatrix( DiscretizationModel.ReservoirGrid );
            
            for ph=1:obj.NofPhases
                %% Transmissibilities of flow and heat
                % Reservoir
                Index.Start = 1;
                Index.End = Nm;
                % Compute (convective) heat transmissibility and gravity
                [obj.Thph{ph, 1}, obj.Ghph{ph, 1}] = obj.MatrixAssembler.ConvectiveHeatTransmissibilityMatrix( ...
                    DiscretizationModel.ReservoirGrid, ...
                    obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, ...
                    ProductionSystem.Reservoir.State.Properties(['h_',num2str(ph)]).Value, ...
                    obj.GravityModel.RhoInt{ph, 1} );
            end
            %% Energy Balance Residual
            % Reservoir
            Index.Start = 1;
            Index.End = Nm;
            R_EB = BuildMediumHeatResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, Qhw, Qhf, RTf, 0);
                            
            %% Mass Balance Jacobian (J_MB_P)
            Nm = DiscretizationModel.ReservoirGrid.N;
            Reservoir = ProductionSystem.Reservoir;
            Wells = ProductionSystem.Wells;
            Index.Start = 1;
            Index.End = Nm;
            % Heat Jacobian blocks
            [~, J_EB_H] = BuildMediumHeatJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0); % 0 = reservoir

            % Solve for deltaP
            deltaH = J_EB_H\(-R_EB);
            
            %% Update Reservoir State
            Hm = ProductionSystem.Reservoir.State.Properties('hTfluid');
            Hm.update(deltaH(1:Nm));
            % Update ALL properties
            obj.ComputeProperties(ProductionSystem, FluidModel)            
            % 7. Update wells
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
        end       
        function delta = UpdateState(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                return
            else
%                 delta = obj.UpdateEnthalpy(delta, ProductionSystem, FluidModel, DiscretizationModel);

                Nm = DiscretizationModel.ReservoirGrid.N;
                Nt = DiscretizationModel.N;
                deltaP = delta(1:Nt);
                deltaH = delta(Nt+1:2*Nt);
                
                %% Update reservoir state
                % 1. Update Pressure
                Pm = ProductionSystem.Reservoir.State.Properties('P_2');
                Pm.update(deltaP(1:Nm));
                % 2. Update Enthalpy
                Hm = ProductionSystem.Reservoir.State.Properties('hTfluid');
                Hm.update(deltaH(1:Nm));
                                
                % Update properties and derivatives
                obj.ComputeProperties(ProductionSystem, FluidModel)            
                % 7. Update wells
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
                        H = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('hTfluid');
                        H.update(deltaH(IP:EP));
                        
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
        function [Qw, Qhw]= ComputeSourceTerms(obj, N, Wells)
            Qw = zeros(N, obj.NofPhases);
            Qhw = zeros(N, obj.NofPhases);  
            %Injectors
            for i=1:Wells.NofInj
                c = Wells.Inj(i).Cells;
                Qw(c, :) = Wells.Inj(i).QPhases(:,:); % MB source term
                Qhw(c, :) = Wells.Inj(i).Qh(:,:);     % EB source term 
                % function is UpdateState in
                % ProductionSystem>Wells>injector_pressure.m [class]
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
        function [J_MB_P, J_MB_H] = AddWellsToMassBalanceJacobian(obj, J_MB_P, J_MB_H, State, Wells, K)
            % Define Local handles
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                dQdp = Inj(i).ComputeWellMassFluxDerivativeWithRespectToPressure(K, obj.NofPhases);
                dQdh = Inj(i).ComputeWellMassFluxDerivativeWithRespectToEnthalpy(obj.NofPhases);
                for ph=1:obj.NofPhases
                    for j=1:length(a)
                        % add derivative of inj well to J_MB_P
                        J_MB_P(a(j),a(j)) = J_MB_P(a(j),a(j)) - dQdp(j, ph);
                        J_MB_H(a(j),a(j)) = J_MB_H(a(j),a(j)) - dQdh(j, ph);
                    end
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                dQdp = Prod(i).ComputeWellMassFluxDerivativeWithRespectToPressure(State, K, obj.Mob, obj.drhodp, obj.dMobdp, obj.NofPhases);
                dQdh = Prod(i).ComputeWellMassFluxDerivativeWithRespectToEnthalpy(State, K, obj.Mob, obj.drhodh, obj.dMobdh, obj.NofPhases);
                for ph=1:obj.NofPhases
                    for j=1:length(b)
                        % add derivative of prod well to J_MB_P & J_MB_H
                        J_MB_P(b(j),b(j)) = J_MB_P(b(j),b(j)) - dQdp(j, ph);
                        J_MB_H(b(j),b(j)) = J_MB_H(b(j),b(j)) - dQdh(j, ph);
                    end
                end
            end
        end
        function [J_EB_P, J_EB_H] = AddWellsToEnergyBalanceJacobian(obj, J_EB_P, J_EB_H, State, Wells, K)
            % Define Local handles
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                dQhdp = Inj(i).ComputeWellHeatFluxDerivativeWithRespectToPressure(K, obj.NofPhases);
                dQhdh = Inj(i).ComputeWellHeatFluxDerivativeWithRespectToEnthalpy(obj.NofPhases);
                for ph=1:obj.NofPhases
                    for j=1:length(a)
                        % add derivative of inj well to J_EB_P
                        J_EB_P(a(j),a(j)) = J_EB_P(a(j),a(j)) - dQhdp(j, ph);
                        J_EB_H(a(j),a(j)) = J_EB_H(a(j),a(j)) - dQhdh(j, ph);
                    end
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                dQhdp = Prod(i).ComputeWellHeatFluxDerivativeWithRespectToPressure(State, K, obj.Mob, obj.drho_times_hdp, obj.dMobdp, obj.NofPhases);
                dQhdh = Prod(i).ComputeWellHeatFluxDerivativeWithRespectToEnthalpy(State, K, obj.Mob, obj.drhodh, obj.dMobdh, obj.NofPhases);
                for ph=1:obj.NofPhases
                    for j=1:length(b)
                        % add derivative of prod well to J_EB_P & J_EB_H
                        J_EB_P(b(j),b(j)) = J_EB_P(b(j),b(j)) - dQhdp(j, ph);
                        J_EB_H(b(j),b(j)) = J_EB_H(b(j),b(j)) - dQhdh(j, ph);
                    end
                end
            end
        end
        function d2Fdh2 = ComputeSecondDerivativeHeatConvectionFlux(obj, Medium, FluidModel, hTfluid)
            d2Fdh2 = zeros(size(obj.Mob,1),1);
            
            % Get P and H indices 
            [~,Pindex] = ismember( round(Medium.State.Properties('P_2').Value,-log10(FluidModel.Pstepsize)), round(FluidModel.Ptable,-log10(FluidModel.Pstepsize)) );
            [~,Hindex] = ismember( round(hTfluid,-log10(FluidModel.Hstepsize)), round(FluidModel.Htable,-log10(FluidModel.Hstepsize)) );

            STable = FluidModel.TablePH.('S_1'); 
            S1 = FluidModel.Phases(1).GetSaturation(Pindex, Hindex, STable); 
            kr = FluidModel.RelPermModel.ComputeRelPerm(FluidModel.Phases, S1);

            for i=1:obj.NofPhases
                % Density
                rhoTable = FluidModel.TablePH.(['rho_', num2str(i)]);
                rho(:,i) = FluidModel.Phases(i).GetDensity(Pindex, Hindex, rhoTable);
                drhodh(:,i) = FluidModel.Phases(i).ComputeDrhoDh(Pindex, Hindex, rhoTable);
                
                % Phase Enthalpy 
                PhaseEnthalpyTable = FluidModel.TablePH.(['H_',num2str(i)]);
                h(:,i) = FluidModel.Phases(i).GetPhaseEnthalpy(Pindex, PhaseEnthalpyTable);    

                % 2nd derivative density wrt enthalpy
                d2rhodh2(:,i) = FluidModel.Phases(i).ComputeD2rhoDh2(Pindex, Hindex, rhoTable);
                
                % Viscosity
                muTable = FluidModel.TablePH.(['mu_',num2str(i)]);
                mu(:,i) = FluidModel.Phases(i).GetViscosity(Pindex, Hindex, muTable);
                
                % Mobility
                Mob(:,i) = kr(:,i)./mu(:,i);

                % 1st derivatives Viscosity and Mobility
                dmudh(:,i) = FluidModel.Phases(i).ComputeDmuDh(Pindex, Hindex, muTable);
                dMobdh(:,i) = ( -1 .* dmudh(:,i) .* kr(:,i) ) ./ mu(:,i).^2;
                
                % 2nd derivatives Viscosity and Mobility
                d2mudh2(:,i) = FluidModel.Phases(i).ComputeD2muDh2(Pindex, Hindex, muTable);
                d2Mobdh2(:,i) = -1 .* kr(:,i) .* ( d2mudh2(:,i) ./ mu(:,i).^2 - 2 .* dmudh(:,i) ./ mu(:,i).^3 );   
            end
            
            % Correction to Mobility and its derivatives
            Mob(isnan(Mob))=0; 
            dMobdh(isnan(dMobdh))=0;
            d2Mobdh2(isnan(d2Mobdh2))=0;

            % F = rho .* h .* Mob;
            % dFdh = drhodh .* h .* Mob + rho .* h .* dMobdh;
            
            % d2Fd2h = d2rhod2h .* h .* Mob + drhodh .* h .* dMobdh + drhodh .* h .* dMobdh + rho .* h .* d2Mobd2h;
            % per phase
            for ph=1:obj.NofPhases
                d2Fdh2 = d2Fdh2 + ...
                         d2rhodh2(:,ph) .* h(:,ph) .* Mob(:,ph) ;%    + ...
%                          drhodh(:,ph)   .* h(:,ph) .* dMobdh(:,ph) + ...
%                          drhodh(:,ph)   .* h(:,ph) .* dMobdh(:,ph) + ...
%                          rho(:,ph)      .* h(:,ph) .* d2Mobdh2(:,ph) ;            
            end 
        end
        function delta = UpdateEnthalpy(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            % Detection of inflexion point (you have a bunch of zeros...)
            Nt = DiscretizationModel.N;
            deltaP_Old = delta(1:Nt);
            deltaH_Old = delta(Nt+1:2*Nt);

            H_old = ProductionSystem.Reservoir.State.Properties('hTfluid').Value;
            H_new = H_old + deltaH_Old;
            
            % Flux Correction - Mousa and Me
            d2Fdh2_old = obj.ComputeSecondDerivativeHeatConvectionFlux( ProductionSystem.Reservoir, FluidModel, H_old);
            d2Fdh2_new = obj.ComputeSecondDerivativeHeatConvectionFlux( ProductionSystem.Reservoir, FluidModel, H_new);
            H_new = H_new.*(d2Fdh2_new.*d2Fdh2_old >= 0) + 0.5*(H_new + H_old).*(d2Fdh2_new.*d2Fdh2_old < 0);
            deltaH_New = H_new - H_old;
            
            % WE ARE HERE, THE PART ABOVE IS IMPLEMENTED
            
            % This is the actual update of the delta solution (guess)
            deltaP_New = deltaP_Old .* (deltaH_New./deltaH_Old); 
            % Not sure about this part (line 823); it appears that the pressure is/must also be corrected,
            % i.e. steered into the same direction as enthalpy (formerly temperature)
            delta(   1:  Nt) = deltaP_New;
            delta(Nt+1:2*Nt) = deltaH_New;
            
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