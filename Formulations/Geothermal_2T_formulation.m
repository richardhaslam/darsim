% Geothermal formulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Geothermal_2T_formulation < formulation
    properties
        drho % first column is drhodp second column is drhodT
        dmu
        dh % first column is dhdp second column is dhdT
        Cp % fluid specific heat
        Kf % fluid conductivity
        AveragedTemperature
    end
    methods
        function obj = Geothermal_2T_formulation()
            obj@formulation();
            obj.Tph = cell(2,1);
            obj.Gph = cell(2,1);
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
                 x(Nt+Index.Start  :Nt + Index.End  ) = ProductionSystem.Reservoir.State.Properties('Tf').Value;
                 x(2*Nt+Index.Start:2*Nt + Index.End) = ProductionSystem.Reservoir.State.Properties('Tr').Value;
                 
                 % Get fracture unknowns
                 for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                     Index.Start = Index.End+1;
                     Index.End = Index.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                     x(Index.Start     :Index.End       ) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_',num2str(ph)]).Value;
                     x(Nt+Index.Start  :Nt + Index.End  ) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('Tf').Value;
                     x(2*Nt+Index.Start:2*Nt + Index.End) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('Tr').Value;   
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
            %% 1. Reservoir Properties Only for Geothermal
            obj.Cp = FluidModel.Phases.Cp; % define Cp from simualtion builder
            obj.Kf = FluidModel.Phases.Kf; % define kf from simualtion builder
            %% 2. Reservoir Properties and Derivatives
            obj.drho = FluidModel.ComputeDrho(ProductionSystem.Reservoir.State); % call derivative of density from geothermal fluid model
            obj.dh = FluidModel.ComputeDh(ProductionSystem.Reservoir.State); % call derivative of enthalpy from geothermal fluid model
            obj.dmu = FluidModel.ComputeDmu(ProductionSystem.Reservoir.State); % call derivative of viscosity from geothermal fluid model
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.Properties('mu_1').Value); 
            obj.dMob = FluidModel.ComputeDMobdT(ProductionSystem.Reservoir.State);
            %% 3. Fractures Properties and Derivatives
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                obj.drho = [obj.drho; FluidModel.ComputeDrho(ProductionSystem.FracturesNetwork.Fractures(f).State) ];
                obj.dh = [obj.dh; FluidModel.ComputeDh(ProductionSystem.FracturesNetwork.Fractures(f).State) ];
                obj.dmu = [obj.dmu; FluidModel.ComputeDmu(ProductionSystem.FracturesNetwork.Fractures(f).State) ];
                obj.Mob = [obj.Mob; FluidModel.ComputePhaseMobilities(ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('mu_1').Value)];
                obj.dMob = [obj.dMob; FluidModel.ComputeDMobdT(ProductionSystem.FracturesNetwork.Fractures(f).State)];
            end
        end
        %% Methods for FIM Coupling
        function Residual_P  = BuildMediumFlowResidual(obj, Medium, Grid, dt, State0, Index, qw, qf, f, ph)
            % Create local variables
            N = Grid.N;
            
            % Copy values in local variables
            s_old = State0.Properties(['S_', num2str(ph)]).Value(Index.Start:Index.End);
            rho_old = State0.Properties(['rho_', num2str(ph)]).Value(Index.Start:Index.End);
            P_old = State0.Properties(['P_', num2str(ph)]).Value(Index.Start:Index.End); % P at previous time step
            P = Medium.State.Properties(['P_', num2str(ph)]).Value;
            s = Medium.State.Properties(['S_', num2str(ph)]).Value;
            rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
            depth = Grid.Depth;
            Medium.ComputePorosity(P_old);
            pv_old = Medium.Por*Grid.Volume;
            Medium.ComputePorosity(P);
            pv = Medium.Por*Grid.Volume;
            
            % Accumulation Term
            AS = 1/dt;

            % RESIDUAL
            Residual_P  = AS.*(pv.*rho - pv_old.*rho_old)...
                + obj.Tph{ph, 1+f} * P...
                - obj.Gph{ph, 1+f} * depth...
                - qw(Index.Start:Index.End, ph)...
                - qf(Index.Start:Index.End, ph);
        end
        function Residual_Tf = BuildMediumFluidHeatResidual(obj, Medium, Grid, dt, State0, Index, qhw, qhf, QTfTr, f, ph)
            % Create local variables
            N = Grid.N;
            
            % Copy values in local variables
            s_old = State0.Properties(['S_', num2str(ph)]).Value(Index.Start:Index.End);
            rho_old = State0.Properties(['rho_', num2str(ph)]).Value(Index.Start:Index.End);
            P_old = State0.Properties(['P_', num2str(ph)]).Value(Index.Start:Index.End); % P at previous time step
            P = Medium.State.Properties(['P_', num2str(ph)]).Value;
            Tf_old = State0.Properties('Tf').Value(Index.Start:Index.End); % T at previous time step
            Tf = Medium.State.Properties('Tf').Value;
            Tr = Medium.State.Properties('Tr').Value;
            s = Medium.State.Properties(['S_', num2str(ph)]).Value;
            rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
            depth = Grid.Depth;
            Medium.ComputePorosity(P_old);
            pv_old = Medium.Por*Grid.Volume;
            Medium.ComputePorosity(P);
            pv = Medium.Por*Grid.Volume;
            
            % Accumulation Term
            AS = obj.Cp/dt;
            U = ComputeOverallHeat(obj, Medium, Grid, Index, ph, f);
            area = ComputeArea(obj, Medium, Grid);
            
            % RESIDUAL
            Residual_Tf  = AS.*(pv.*rho.*Tf - pv_old.*rho_old.*Tf_old)...
                + obj.Th{ph, 1+f} * P...
                - U .* area .*(Tr-Tf)...
                - obj.Gph{ph, 1+f} * depth...
                - qhw(Index.Start:Index.End, ph)...
                - qhf(Index.Start:Index.End, ph)...
                - QTfTr(Index.Start:Index.End, ph);
           
        end
        function Residual_Tr = BuildMediumRockHeatResidual(obj, Medium, Grid, dt, State0, Index, QTfTr, f, ph)
            % Create local variables
            N = Grid.N;
            
            % Copy values in local variables
            P_old = State0.Properties(['P_', num2str(ph)]).Value(Index.Start:Index.End); % P at previous time step
            P = Medium.State.Properties(['P_', num2str(ph)]).Value;
            Tr_old = State0.Properties('Tr').Value(Index.Start:Index.End); % T at previous time step
            Tr = Medium.State.Properties('Tr').Value;
            T = Medium.State.Properties('Tf').Value;
            depth = Grid.Depth;
            Medium.ComputePorosity(P_old);
            mv_old = (1 - Medium.Por) * Grid.Volume; % matrix volume
            Medium.ComputePorosity(P);
            mv = (1 - Medium.Por) * Grid.Volume; % matrix volume
            
            % Accumulation Term
            AS = (Medium.Cpr * Medium.Rho)/dt;
            
            U = ComputeOverallHeat(obj, Medium, Grid, Index, ph, f);
            area = ComputeArea(obj, Medium, Grid);
            
            % RESIDUAL
            Residual_Tr  = AS.*(mv.*Tr - mv_old.*Tr_old)...
                + obj.Tk * Tr + U .* area .*(Tr-T)... % Tk = transmisibility of rock conductivity
                - QTfTr(Index.Start:Index.End);
        end
        function ResidualFull= BuildResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Compute source terms
            [Qw, Qhw] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            Qf = zeros(DiscretizationModel.N, obj.NofPhases);
            Qhf= zeros(DiscretizationModel.N, obj.NofPhases);
            QTfTr = zeros(DiscretizationModel.N, obj.NofPhases);
            if ProductionSystem.FracturesNetwork.Active
                [Qf, Qhf, QTfTr] = obj.ComputeSourceTerms_frac_mat(ProductionSystem, DiscretizationModel);
            end

            % Initialise residual vector (Nph * N, 1)
            Nm = DiscretizationModel.ReservoirGrid.N;
            if ProductionSystem.FracturesNetwork.Active
                Nf = DiscretizationModel.FracturesGrid.N;
            else
                Nf = 0;
            end
            Nt = DiscretizationModel.N;
            Residual = zeros( 3*Nm + 2*sum(Nf) , 1 );
                        
            for ph=1:obj.NofPhases
                %% Flow Residual of reservoir
                % Transmissibility of flow for reservoir
                [obj.Tph{ph, 1}, obj.Gph{ph, 1}, ~] = ...
                    obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, obj.GravityModel.RhoInt{ph, 1});
                % Transmissibility of heat for reservoir
                [obj.Th{ph, 1}, ~, ~] = ...
                        obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, obj.UpWind{ph, 1}, obj.Mob(1:Nm, ph), ...
                        ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value .* ...
                        ProductionSystem.Reservoir.State.Properties(['h_',num2str(ph)]).Value, obj.GravityModel.RhoInt{ph, 1});
                Index.Start = 1;
                Index.End = Nm;
                Residualm = BuildMediumFlowResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, Qw, Qf, 0, ph);
                Residual((ph-1)*Nt + Index.Start: (ph-1)*Nt + Index.End) = Residualm;
                %% Flow Residual of fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End+1;
                    Index.End = Index.Start + Nf(f) - 1;
                    % Transmissibility of flow for fractures
                    [obj.Tph{ph, 1+f}, obj.Gph{ph, 1+f}, ~] = ...
                        obj.TransmissibilityMatrix(DiscretizationModel.FracturesGrid.Grids(f), obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph),...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value, obj.GravityModel.RhoInt{ph, 1+f});
                     % Transmissibility of heat for fractures
                    [obj.Th{ph, 1+f}, ~, ~] = ...
                        obj.TransmissibilityMatrix(DiscretizationModel.FracturesGrid.Grids(f), obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph),...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value .* ...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['h_',num2str(ph)]).Value, obj.GravityModel.RhoInt{ph, 1+f});
                    Residualf = BuildMediumFlowResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, Qw, Qf, f, ph);
                    Residual(Index.Start:Index.End) = Residualf;
                end
                
                %% Fluid Heat Residual of reservoir
                ph = 1;
                Index.Start = Index.End + 1;
                Index.End = Index.Start + Nm - 1;
                Index_r.Start = 1;
                Index_r.End = Nm;
                Residualm = BuildMediumFluidHeatResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index_r, Qhw, Qhf, 0*QTfTr, 0, ph);
                Residual(Index.Start: Index.End) = Residualm;
                %% Fluid Heat Residual of fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End + 1;
                    Index.End = Index.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                    Index_r.Start = Index_r.End + 1;
                    Index_r.End = Index_r.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                    Residualf = BuildMediumFluidHeatResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index_r, Qhw, Qhf, QTfTr, f, ph);
                    Residual(Index.Start:Index.End) = Residualf; 
                end
            end
            
            %% Rock Heat Residual of reservoir
            % Transmissibility of rock heat for reservoir
            [~, ~, obj.Tk] = ...
                obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, obj.UpWind{ph, 1}, obj.Mob(1:DiscretizationModel.ReservoirGrid.N, ph), ...
                ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, obj.GravityModel.RhoInt{ph, 1});
            Index_r.Start = 1;
            Index_r.End = Nm;
            Residualm = BuildMediumRockHeatResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index_r, QTfTr, 0, ph);
            Index.Start = Index.End + 1;
            Index.End = Index.Start + Nm - 1;
            Residual(Index.Start:Index.End) = Residualm;
            
            % Option to sum-up two equations of "Tf" and "Tr" and
            % solving for one average temperature.
            if obj.AveragedTemperature == "Off"
                ResidualFull = Residual;
            else % Temperature Tf and Tr will be averaged and two equations will be summed
                ResidualFull = zeros( 2*Nm + 2*sum(Nf) , 1 );
                ResidualFull(1:2*Nt) = Residual( 1:2*Nt );
                ResidualFull(Nt+1:Nt+Nm) = ResidualFull(Nt+1:Nt+Nm) + Residual(2*Nt+1:end);
            end
        end
        function U = ComputeOverallHeat(obj, Medium, Grid, Index, ph, f)
            Ks = mean(Medium.k_cond,2);
            Nu = ComputeNusselt(obj, Medium, Grid, Index, ph, f);
            U = 1./(Medium.Dp./(Nu*obj.Kf) + Medium.Dp./(10*Ks));
            function Nu = ComputeNusselt(obj, Medium, Grid, Index, ph, f)
                Por = Medium.Por;
                rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
                Pr = obj.Cp./(obj.Mob(Index.Start:Index.End) * obj.Kf);
                v = ComputeVelocity(obj, Medium, Index, Grid, f);
                if f ~= 0
                    Re = rho .* v(:) .* Medium.Dp .* obj.Mob(Index.Start:Index.End);
                else
                    Re = 1./(1-Por) .* rho .* v(:) .* Medium.Dp .* obj.Mob(Index.Start:Index.End);
                end
                if ~isreal(Re)
                    error('Reynold is not real!');
                end
                if ~isreal(Pr)
                    error('Prandtl is not real!');
                end
                Nu = 0.225./Por .* Pr.^0.33 .* Re.^0.67;
                
                function v = ComputeVelocity(obj, Medium, Index, Grid, f)
                    % Check upwind scheme
                    Mob = reshape(obj.Mob(Index.Start:Index.End), Grid.Nx, Grid.Ny, Grid.Nz);
                    v_x = ( obj.U{1, f+1}.x(1:end-1,:,:) + obj.U{1, f+1}.x(2:end,:,:) ) .* Mob / 2;
                    v_y = ( obj.U{1, f+1}.y(:,1:end-1,:) + obj.U{1, f+1}.y(:,2:end,:) ) .* Mob / 2;
                    v_z = ( obj.U{1, f+1}.z(:,:,1:end-1) + obj.U{1, f+1}.z(:,:,2:end) ) .* Mob / 2;
                    v = (abs(v_x) + abs(v_y) + abs(v_z)) / 3;
                end
            end
        end
        function area = ComputeArea(obj, Medium, Grid)
            area = 6 * (1 - Medium.Por)/Medium.Dp * Grid.Volume;
%             area = 1e-7*area;
        end
        function [J_P_P , J_P_Tf , J_P_Tr ] = BuildMediumFlowJacobian(obj, Medium, Wells, Grid, dt, Index, f, ph)
            % Create local variables
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            
            P = Medium.State.Properties(['P_', num2str(ph)]).Value;
            Tf = Medium.State.Properties(['Tf']).Value;
            Tr = Medium.State.Properties(['Tr']).Value;
            rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
            h = Medium.State.Properties(['h_', num2str(ph)]).Value;
            mu = Medium.State.Properties(['mu_', num2str(ph)]).Value;
            s = Medium.State.Properties(['S_', num2str(ph)]).Value; % why this one always be zero?

            Medium.ComputeDerPorosity(P);
            por = Medium.Por;
            dpor = Medium.DPor;
            
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            % 1.a Pressure Block
            J_P_P = obj.Tph{ph,1+f};
            drhodp = obj.drho(Index.Start:Index.End, 1);
            
            % 1.b: compressibility part
            dMupx = obj.UpWind{ph,1+f}.x*(obj.Mob(Index.Start:Index.End, ph) .* drhodp);
            dMupy = obj.UpWind{ph,1+f}.y*(obj.Mob(Index.Start:Index.End, ph) .* drhodp);
            dMupz = obj.UpWind{ph,1+f}.z*(obj.Mob(Index.Start:Index.End, ph) .* drhodp);
            
            vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
            vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
            vecY1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
            vecY2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
            vecZ1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
            vecZ2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
            acc = Grid.Volume/dt .* (por .* drhodp + rho .*dpor); 
            
            DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J_P_P = J_P_P + spdiags(DiagVecs, DiagIndx, N, N);
            
            % 2. Fluid Temperature Block
            drhodT = obj.drho(Index.Start:Index.End, 2);
            dMupx = obj.UpWind{ph,1+f}.x * (obj.dMob(Index.Start:Index.End, ph) .* rho + obj.Mob(Index.Start:Index.End, ph) .* drhodT);
            dMupy = obj.UpWind{ph,1+f}.y * (obj.dMob(Index.Start:Index.End, ph) .* rho + obj.Mob(Index.Start:Index.End, ph) .* drhodT);
            dMupz = obj.UpWind{ph,1+f}.z * (obj.dMob(Index.Start:Index.End, ph) .* rho + obj.Mob(Index.Start:Index.End, ph) .* drhodT);
            % Construct JPT block
            x1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
            x2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
            y1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
            y2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
            z1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
            z2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
            v =  Grid.Volume/dt .* por .* drhodT ;
            DiagVecs = [-z2, -y2, -x2, z2+y2+x2-z1-y1-x1+v, x1, y1, z1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J_P_Tf = spdiags(DiagVecs,DiagIndx,N,N);
            
            %Add capillarity
%             JT = JT - Jp * spdiags(obj.dPc, 0, N, N);
            
            % Add Wells
            % for now, we will consider an only 2-phase system for adding the wells to the jacobian
            if f == 0 % only for reservoir
                [J_P_P, J_P_Tf, ~, ~] = obj.AddWellsToJacobian(J_P_P, J_P_Tf, Medium.State, Wells, Medium.K(:,1), ph);
            end
            
            % Jacobian J_pTr
            J_P_Tr = 0*J_P_Tf;
            
        end
        function [J_Tf_P, J_Tf_Tf, J_Tf_Tr] = BuildMediumFluidHeatJacobian(obj, Medium, Wells, Grid, dt, Index, f, ph)
            % Create local variables
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            
            P = Medium.State.Properties(['P_', num2str(ph)]).Value;
            T = Medium.State.Properties(['Tf']).Value;
            Tr = Medium.State.Properties(['Tr']).Value;
            rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
            h = Medium.State.Properties(['h_', num2str(ph)]).Value;
            mu = Medium.State.Properties(['mu_', num2str(ph)]).Value;
            s = Medium.State.Properties(['S_', num2str(ph)]).Value; 

            Medium.ComputeDerPorosity(P);
            por = Medium.Por;
            dpor = Medium.DPor;
            
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            % 1.a Pressure Block
            J_Tf_P = obj.Th{ph, 1+f};
            drhodp = obj.drho(Index.Start:Index.End, 1);
            dhdp = obj.dh(Index.Start:Index.End, 1);
            
            % 1.b: compressibility part
            dMupx = obj.UpWind{ph,1+f}.x*(obj.Mob(Index.Start:Index.End, ph) .* ( drhodp .* h + dhdp .* rho ));
            dMupy = obj.UpWind{ph,1+f}.y*(obj.Mob(Index.Start:Index.End, ph) .* ( drhodp .* h + dhdp .* rho ));
            dMupz = obj.UpWind{ph,1+f}.z*(obj.Mob(Index.Start:Index.End, ph) .* ( drhodp .* h + dhdp .* rho ));
            
            vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
            vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
            vecY1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
            vecY2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
            vecZ1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
            vecZ2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
            acc = Grid.Volume/dt .* obj.Cp .* T .*(por .* drhodp + rho .*dpor); 
            
            DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J_Tf_P = J_Tf_P + spdiags(DiagVecs, DiagIndx, N, N);
            
            % 2. Fluid Temperature Block
            drhodT = obj.drho(Index.Start:Index.End, 2);
            dhdT = obj.dh(Index.Start:Index.End, 2);
            dMob = obj.dMob(Index.Start:Index.End, ph);
            Mob = obj.Mob(Index.Start:Index.End, ph);
            U = ComputeOverallHeat(obj, Medium, Grid, Index, ph, f);
            area = ComputeArea(obj, Medium, Grid);
            
            dMupx = obj.UpWind{ph,1+f}.x * (dMob .* rho .* h + drhodT .* Mob .* h  + dhdT .* rho .* Mob);
            dMupy = obj.UpWind{ph,1+f}.y * (dMob .* rho .* h + drhodT .* Mob .* h  + dhdT .* rho .* Mob);
            dMupz = obj.UpWind{ph,1+f}.z * (dMob .* rho .* h + drhodT .* Mob .* h  + dhdT .* rho .* Mob);
            % Construct JPT block
            x1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0)   .* dMupx;
            x2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
            y1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0)   .* dMupy;
            y2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
            z1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0)   .* dMupz;
            z2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
            v1 = Grid.Volume/dt .* obj.Cp .* por .*( drhodT .* T + rho );
            v2 = U .* area;
            DiagVecs = [-z2, -y2, -x2, z2+y2+x2-z1-y1-x1+v1+v2, x1, y1, z1];
%             DiagVecs = [-z2, -y2, -x2, z2+y2+x2-z1-y1-x1+v1, x1, y1, z1]; 
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J_Tf_Tf = spdiags(DiagVecs,DiagIndx,N,N);
            
            % Add Wells
            % for now, we will consider an only 2-phase system for adding the wells to the jacobian
            if f == 0 % only for reservoir
                [~, ~, J_Tf_P, J_Tf_Tf] = obj.AddWellsToJacobian(J_Tf_P, J_Tf_Tf, Medium.State, Wells, Medium.K(:,1), ph);
            end
            
            % 3. Rock Temperature Block
            v3 = -U .* area;
            DiagVecs = [-z2*0, -y2*0, -x2*0, v3, x1*0, y1*0, z1*0];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J_Tf_Tr = spdiags(DiagVecs,DiagIndx,N,N);
            
        end
        function [J_Tr_P, J_Tr_Tf, J_Tr_Tr] = BuildMediumRockHeatJacobian(obj, Medium, Grid, dt, Index, f, ph)
            % Create local variables
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            
            P = Medium.State.Properties(['P_', num2str(ph)]).Value;
            T = Medium.State.Properties(['Tf']).Value;
            Tr = Medium.State.Properties(['Tr']).Value;
            rho = Medium.State.Properties(['rho_', num2str(ph)]).Value;
            h = Medium.State.Properties(['h_', num2str(ph)]).Value;
            mu = Medium.State.Properties(['mu_', num2str(ph)]).Value;
            s = Medium.State.Properties(['S_', num2str(ph)]).Value; 

            Medium.ComputeDerPorosity(P);
            por = Medium.Por;
            dpor = Medium.DPor;
            
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            % 1. Pressure Block
            drhodp = obj.drho(Index.Start:Index.End, 1);
            dhdp = obj.dh(Index.Start:Index.End, 1);
            
            vecX1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0);
            vecX2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0);
            vecY1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0);
            vecY2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0);
            vecZ1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0);
            vecZ2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0);
            acc = Grid.Volume/dt .* Medium.Cpr .* Medium.Rho .* Tr .* -dpor; 
            
            DiagVecs = [-vecZ2*0, -vecY2*0, -vecX2*0, acc, vecX1*0, vecY1*0, vecZ1*0];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J_Tr_P = spdiags(DiagVecs, DiagIndx, N, N);
            
            % 2. Fluid Temperature Block
            drhodT = obj.drho(Index.Start:Index.End, 2);
            dhdT = obj.dh(Index.Start:Index.End, 2);
            dMob = obj.dMob(Index.Start:Index.End, ph);
            Mob = obj.Mob(Index.Start:Index.End, ph);
            U = ComputeOverallHeat(obj, Medium, Grid, Index, ph, f);
            area = ComputeArea(obj, Medium, Grid);
            
            % Construct JPT block
            x1 = min(reshape(obj.U{ph,1+f}.x(1:Nx,:,:), N, 1), 0);
            x2 = max(reshape(obj.U{ph,1+f}.x(2:Nx+1,:,:), N, 1), 0);
            y1 = min(reshape(obj.U{ph,1+f}.y(:,1:Ny,:), N, 1), 0);
            y2 = max(reshape(obj.U{ph,1+f}.y(:,2:Ny+1,:), N, 1), 0);
            z1 = min(reshape(obj.U{ph,1+f}.z(:,:,1:Nz), N, 1), 0);
            z2 = max(reshape(obj.U{ph,1+f}.z(:,:,2:Nz+1), N, 1), 0);
            v1 = -U .* area;
            DiagVecs = [-z2*0, -y2*0, -x2*0, v1, x1*0, y1*0, z1*0];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J_Tr_Tf = spdiags(DiagVecs,DiagIndx,N,N);
            
            % 3. Rock Temperature Block
            J_Tr_Tr = obj.Tk;
            v2 = Grid.Volume/dt .* (1 - por) .* Medium.Cpr .* Medium.Rho;
            v3 = U .* area;
            DiagVecs = [-z2*0, -y2*0, -x2*0, v2+v3, x1*0, y1*0, z1*0];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J_Tr_Tr = J_Tr_Tr + spdiags(DiagVecs,DiagIndx,N,N);
            
        end
        function JacobianFull = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            %% Jacobian's assembly for single phase geothermal reservoir
            % | JPP  JPTf  JPTr  | dP  |    | RP  | <-- 3 & 2 unknowns --> | JPP  JPT | dP |    | RP |
            % | JTfP JTfTf JTfTr | dTf | = -| RTf |                        | JTP  JTT | dT | = -| RT |
            % | JTrP JTrTf JTrTr | dTr |    | RTr |
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
                [J_P_P, J_P_Tf, J_P_Tr] = BuildMediumFlowJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0, ph); % 0 = reservoir , 1 = phase
                
                % Heat Jacobian blocks
                [J_Tf_P, J_Tf_Tf, J_Tf_Tr] = BuildMediumFluidHeatJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0, ph); % 0 = reservoir , 1 = phase
                
                % Rock Heat Jacobian blocks
                [J_Tr_P, J_Tr_Tf, J_Tr_Tr] = BuildMediumRockHeatJacobian(obj, Reservoir, DiscretizationModel.ReservoirGrid, dt, Index, 0, 1); % 0 = reservoir , 1 = phase
                
                %% Jacobian of the fractures
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Nf = DiscretizationModel.FracturesGrid.N;
                    Index.Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                    Index.End = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, Nf(f));
                    
                    % Flow Jacobian Block
                    [Jf_PP, Jf_PT, ~] = BuildMediumFlowJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f, ph);
                    J_P_P  = blkdiag(J_P_P, Jf_PP);
                    J_P_Tf  = blkdiag(J_P_Tf, Jf_PT);
                    %J_PTr  = blkdiag(J_PTr, Jf_PTr); fractures do not have Tr
                    J_P_Tr = vertcat(J_P_Tr, sparse(Nf(f), Nm));
                    
                    % Heat Jacobian Block
                    [Jf_TP, Jf_TT, ~] = BuildMediumFluidHeatJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index, f, ph);
                    J_Tf_P  = blkdiag(J_Tf_P, Jf_TP);
                    J_Tf_Tf  = blkdiag(J_Tf_Tf, Jf_TT);
                    %J_TTr  = blkdiag(J_TTr, Jf_TTr); fractures do not have Tr
                    J_Tf_Tr = vertcat(J_Tf_Tr, sparse(Nf(f), Nm));
                    
                    % Rock Heat Jacobian Block
                    %J_TrP  = blkdiag(J_TrP, 0*Jf_TP);   fractures do not have Tr
                    J_Tr_P = horzcat(J_Tr_P, sparse(Nm, Nf(f)));
                    %J_TrT  = blkdiag(J_TrT, 0*Jf_TT);   fractures do not have Tr
                    J_Tr_Tf = horzcat(J_Tr_Tf, sparse(Nm, Nf(f)));
                    %J_TrTr = blkdiag(J_TrTr, 0*Jf_TTr); fractures do not have Tr
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
                k_cond = ProductionSystem.Reservoir.k_cond(:,1);
                obj.drhodp = obj.drho(:,1);
                drhodT = obj.drho(:,2);
                dhdp = obj.dh(:,1);
                dhdT = obj.dh(:,2);
                for c = 1:length(DiscretizationModel.CrossConnections)
                    T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                    UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                    i = c + Nm;
                    j = DiscretizationModel.CrossConnections(c).Cells;
                    jj = j(j<=Nm); % a temporary "j" containing of only matrix indeces;
                    
                    %% J_P_P Coupling
                    J_P_P_1_conn = - T_Geo .*  UpWind(:, ph) .* obj.Mob(j, ph) .* ( rho(j,ph) + obj.drhodp(j, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_P_P_2_conn = - T_Geo .* ~UpWind(:, ph) .* obj.Mob(i, ph) .* ( rho(i,ph) + obj.drhodp(i, ph) .* ( P(j,ph) - P(i,ph) ) );
                    J_P_P_conn  = J_P_P_1_conn + J_P_P_2_conn;
                    
                    % frac - mat or frac1 - frac2
                    if sum(J_P_P(i, j))~=0 || sum(J_P_P(j, i))~=0
                        error('J_P_P(i, j) or J_P_P(j, i) is not zero!');
                    end
                    J_P_P(i, j) = J_P_P_conn;
                    J_P_P(i, i) = J_P_P(i, i) - sum(J_P_P_conn);
                    % mat-frac or frac2 - frac1
                    J_P_P(j, i) = J_P_P_conn';
                    J_P_P(sub2ind([Nt, Nt], j, j)) = J_P_P(sub2ind([Nt, Nt], j, j)) - J_P_P_conn;
                    
                    %% J_P_Tf Coupling
                    J_P_Tf_1_conn = - T_Geo .*  UpWind(:, ph) .* (P(j, ph) - P(i, ph)).* ( rho(j, ph) .* obj.dMob(j, ph) + obj.Mob(j, ph) .* drhodT(j, ph) );
                    J_P_Tf_2_conn = - T_Geo .* ~UpWind(:, ph) .* (P(j, ph) - P(i, ph)).* ( rho(i, ph) .* obj.dMob(i, ph) + obj.Mob(i, ph) .* drhodT(i, ph) );
                    
                    % frac - mat or frac1 - frac2
                    if sum(J_P_Tf(i, j))~=0 || sum(J_P_Tf(j, i))~=0
                        error('J_P_Tf(i, j) or J_P_Tf(j, i) is not zero!');
                    end
                    J_P_Tf(i, j) = J_P_Tf_1_conn;
                    J_P_Tf(i, i) = J_P_Tf(i, i) + sum(J_P_Tf_2_conn);
                    % mat-frac or frac2 - frac1
                    J_P_Tf(j, i) = - J_P_Tf_2_conn;
                    % diag of mat or frac2
                    J_P_Tf(sub2ind([Nt, Nt], j, j)) = J_P_Tf(sub2ind([Nt, Nt], j, j)) - J_P_Tf_1_conn;
                    
                    %% J_P_Tr Coupling
                    % J_P_Tr is zero.
                    
                    %% J_Tf_P Coupling
                    J_Tf_P_1_conn = - T_Geo .*  UpWind(:, ph) .* obj.Mob(j, ph) .* ( (P(j, ph) - P(i, ph)) .* (rho(j, ph) .* dhdp(j, ph) + h(j, ph) .* obj.drhodp(j, ph)) ...
                        + rho(j, ph) .* h(j, ph) );
                    J_Tf_P_2_conn = - T_Geo .* ~UpWind(:, ph) .* obj.Mob(i, ph) .* ( (P(j, ph) - P(i, ph)) .* (rho(i, ph) .* dhdp(i, ph) + h(i, ph) .* obj.drhodp(i, ph)) ...
                        + rho(i, ph) .* h(i, ph) );
                    J_Tf_P_conn = J_Tf_P_1_conn + J_Tf_P_2_conn;
                    % frac - mat or frac1 - frac2
                    if sum(J_Tf_P(i, j))~=0 || sum(J_Tf_P(j, i))~=0
                        error('J_Tf_P(i, j) or J_Tf_P(j, i) is not zero!');
                    end
                    J_Tf_P(i, j) = J_Tf_P_conn;
                    J_Tf_P(i, i) = J_Tf_P(i, i) - sum(J_Tf_P_conn);
                    % mat-frac or frac2 - frac1
                    J_Tf_P(j, i) = J_Tf_P_conn';
                    % diag of mat or frac2
                    J_Tf_P(sub2ind([Nt, Nt], j, j)) = J_Tf_P(sub2ind([Nt, Nt], j, j)) - J_Tf_P_conn;
                    
                    %% J_Tf_Tf Coupling
                    J_Tf_Tf_1_conn = - T_Geo .*  UpWind(:, ph) .* ( (P(j, ph) - P(i, ph)) .* ( obj.dMob(j, ph) .* rho(j, ph) .* h(j, ph) + obj.Mob(j, ph) .* drhodT(j, ph) .* h(j, ph) ...
                        + obj.Mob(j, ph) .* rho(j, ph) .* dhdT(j, ph) ) );
                    J_Tf_Tf_2_conn = - T_Geo .* ~UpWind(:, ph) .* ( (P(j, ph) - P(i, ph)) .* ( obj.dMob(i, ph) .* rho(i, ph) .* h(i, ph) + obj.Mob(i, ph) .* drhodT(i, ph) .* h(i, ph) ...
                        + obj.Mob(i, ph) .* rho(i, ph) .* dhdT(i, ph) ) );
                    % frac - mat or frac1 - frac2
                    if sum(J_Tf_Tf(i, j))~=0 || sum(J_Tf_Tf(j, i))~=0
                        error('J_Tf_Tf(i, j) or J_Tf_Tf(j, i) is not zero!');
                    end
                    J_Tf_Tf(i, j) = J_Tf_Tf_1_conn';
                    J_Tf_Tf(i, i) = J_Tf_Tf(i, i) + sum(J_Tf_Tf_2_conn);
                    % mat-frac or frac2 - frac1
                    J_Tf_Tf(j, i) = -J_Tf_Tf_2_conn;
                    % diag of mat or frac2
                    J_Tf_Tf(sub2ind([Nt, Nt], j, j)) = J_Tf_Tf(sub2ind([Nt, Nt], j, j)) - J_Tf_Tf_1_conn;
                    
                    % frac - mat (Conduction part)
                    J_Tf_Tf_conn = T_Geo(1:length(jj)) .* k_cond(jj);
                    J_Tf_Tf(i, i ) = J_Tf_Tf(i, i) + sum(J_Tf_Tf_conn);
                    % Note that "AU(Tr-Tf)" is not neglected here. It is
                    % taken into account in term " k_cond(jj) .* (-1) "
                    
                    %% J_Tf_Tr Coupling
                    J_Tf_Tr_conn = - T_Geo(1:length(jj)) .* k_cond(jj);
                    if sum(J_Tf_Tr(i, jj))~=0
                        error('J_Tf_Tr(i, jj) is not zero!');
                    end
                    J_Tf_Tr(i , jj) = J_Tf_Tr_conn;
                    %J_Tf_Tr(sub2ind([Nt, Nm], jj, jj)) = J_Tf_Tr(sub2ind([Nt, Nm], jj, jj)) - J_Tf_Tr_conn;
                    
                    %% J_Tr_P
                    % J_Tr_P is zero.
                    
                    %% J_Tr_Tf Coupling
                    J_Tr_Tf_conn = - T_Geo(1:length(jj)) .* k_cond(jj);
                    if sum(J_Tr_Tf(jj, i))~=0
                        error('J_Tf_Tr(jj, i) is not zero!');
                    end
                    J_Tr_Tf (jj , i) = J_Tr_Tf_conn';
                    
                    %% J_Tr_Tr Coupling
                    J_Tr_Tr_conn = T_Geo(1:length(jj)) .* k_cond(jj);
                    J_Tr_Tr(sub2ind([Nm, Nm], jj, jj)) = J_Tr_Tr(sub2ind([Nm, Nm], jj, jj)) + J_Tr_Tr_conn;
                    % The coupling of J_Tr_Tr is non-existent for fractures
                    
                end
            end
            
            % Build & Stack Jacobian
            % Option to sum-up two equations of "Tf" and "Tr" and
            % solving for one average temperature.
            if obj.AveragedTemperature == "Off"
                JacobianFull = [J_P_P, J_P_Tf, J_P_Tr ; J_Tf_P, J_Tf_Tf, J_Tf_Tr ; J_Tr_P, J_Tr_Tf, J_Tr_Tr];
            else % Temperature Tf and Tr will be averaged and two equations will be summed
                J_P_Tf (:,1:Nm)    = J_P_Tf (:,1:Nm)    + J_P_Tr;
                J_Tf_Tf(:,1:Nm)    = J_Tf_Tf(:,1:Nm)    + J_Tf_Tr;
                J_Tf_P (1:Nm,:)    = J_Tf_P (1:Nm,:)    + J_Tr_P;
                J_Tf_Tf(1:Nm,:)    = J_Tf_Tf(1:Nm,:)    + J_Tr_Tf; 
                J_Tf_Tf(1:Nm,1:Nm) = J_Tf_Tf(1:Nm,1:Nm) + J_Tr_Tr;
                JacobianFull = [J_P_P, J_P_Tf ; J_Tf_P, J_Tf_Tf];
            end 
        end
        function ConstrainedPressureResidual(obj, FluidModel, ProductionSystem, DiscretizationModel, dt, State0)
            % Initialise residual vector (Nph * N, 1)
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            Nt = DiscretizationModel.N;
            Nm = DiscretizationModel.ReservoirGrid.N;
            
            Reservoir = ProductionSystem.Reservoir;
            Fractures = ProductionSystem.FracturesNetwork.Fractures;
            Wells = ProductionSystem.Wells;
            Index.Start = 1;
            Index.End = Nm;
            
            %% Build Residual & Jacobian Reservoir
            % Compute source terms
            [Qw, ~] = obj.ComputeSourceTerms(DiscretizationModel.N, ProductionSystem.Wells);
            qf = zeros(DiscretizationModel.N, obj.NofPhases);
            
            for ph=1:obj.NofPhases
                % Transmissibility of flow for reservoir
                [obj.Tph{ph, 1}, obj.Gph{ph, 1}, ~] = ...
                    obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, obj.UpWind{ph, 1}, obj.Mob(1:DiscretizationModel.ReservoirGrid.N, ph), ...
                    ProductionSystem.Reservoir.State.Properties(['rho_',num2str(ph)]).Value, obj.GravityModel.RhoInt{ph, 1});
                
                % Building Reservoir Flow Residual
                Rp = BuildMediumFlowResidual(obj, ProductionSystem.Reservoir, DiscretizationModel.ReservoirGrid, dt, State0, Index, Qw, qf, 0, ph);
                
                % Building Reservoir Flow Jacobian
                [J_PP, ~, ~] = BuildMediumFlowJacobian(obj, Reservoir, Wells, DiscretizationModel.ReservoirGrid, dt, Index, 0, 1); % 0 = reservoir , 1 = phase
                
                %% Build Residual & Jacobian Fracture
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End+1;
                    Index.End = Index.Start + DiscretizationModel.FracturesGrid.N(f) - 1;
                    % Transmissibility of flow for fractures
                    [obj.Tph{ph, 1+f}, obj.Gph{ph, 1+f}, ~] = ...
                        obj.TransmissibilityMatrix(DiscretizationModel.FracturesGrid.Grids(f), obj.UpWind{ph, 1+f}, obj.Mob(Index.Start:Index.End, ph),...
                        ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['rho_',num2str(ph)]).Value, obj.GravityModel.RhoInt{ph, 1+f});
                    % BuildResidual of flow for Fractures
                    Rp_f = BuildMediumFlowResidual(obj, ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), dt, State0, Index, Qw, qf, f, ph);
                    Rp = [Rp ; Rp_f];
                    
                    % Start building Flow Jacobian Fracture
                    Nf = DiscretizationModel.FracturesGrid.N;
                    Index_Fracture.Start = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, 1);
                    Index_Fracture.End = DiscretizationModel.Index_Local_to_Global(Nx, Ny, Nz, f, DiscretizationModel.FracturesGrid.Grids(f).N);
                    
                    % Flow Jacobian Block
                    [Jf_PP, ~, ~] = BuildMediumFlowJacobian(obj, Fractures(f), Wells, DiscretizationModel.FracturesGrid.Grids(f), dt, Index_Fracture, f, ph);
                    J_PP  = blkdiag(J_PP, Jf_PP);
                end
                
                %% ADD frac-matrix and frac-frac connections
                % Global variables
                if ProductionSystem.FracturesNetwork.Active
                    FineGrid = [DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids];
                else
                    FineGrid = DiscretizationModel.ReservoirGrid;
                end
                P = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'P_'); % useful for cross connections assembly
                rho = ProductionSystem.CreateGlobalVariables(FineGrid, obj.NofPhases, 'rho_'); % useful for cross connections assembly
                obj.drhodp = obj.drho(:,1);
                for c = 1:length(DiscretizationModel.CrossConnections)
                    T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                    UpWind = DiscretizationModel.CrossConnections(c).UpWind;
                    i = c + Nm;
                    j = DiscretizationModel.CrossConnections(c).Cells;
                    Jp_conn = - T_Geo .* (...
                             UpWind(:, ph) .* obj.Mob(j, ph) .* (rho(j,ph) + obj.drhodp(j, ph) .* ( P(j,ph) - P(i,ph) )) + ...
                            ~UpWind(:, ph) .* obj.Mob(i, ph) .* (rho(i,ph) + obj.drhodp(i, ph) .* ( P(j,ph) - P(i,ph) )) );
                    % frac - mat or frac1 - frac2
                    J_PP(i, j) = Jp_conn;
                    J_PP(i, i) = J_PP(i, i) - sum(Jp_conn);
                    % mat-frac or frac2 - frac1
                    J_PP(j, i) = Jp_conn';
                    J_PP(sub2ind([Nt, Nt], j, j)) = J_PP(sub2ind([Nt, Nt], j, j)) - Jp_conn;
                end
                
            end
            
            % Solve for deltaP
            deltaP = J_PP\(-Rp);
            Pm = ProductionSystem.Reservoir.State.Properties('P_1');
            Pm.update(deltaP(1:Nm));
            
            %% Update Reservoir State
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
        function delta = UpdateState(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
             if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                return
            else

                Nm = DiscretizationModel.ReservoirGrid.N;
                Nt = DiscretizationModel.N;
                deltaP  = delta(1:Nt);
                deltaTf = delta(Nt+1:2*Nt);
                
                % Option to sum-up two equations of "Tf" and "Tr" and
                % solving for one average temperature.
                if obj.AveragedTemperature == "Off"
                    deltaTr = delta(2*Nt+1:end);
                else % Temperature Tf and Tr will be averaged and two equations will be summed
                    deltaTr = deltaTf(1:Nm,1);
                end
                
                %% Update reservoir state
                % 1. Update Pressure
                Pm = ProductionSystem.Reservoir.State.Properties(['P_1']);
                Pm.update(deltaP(1:Nm));
                % 2. Update Fluid Temperature
                Tfm = ProductionSystem.Reservoir.State.Properties('Tf');
                Tfm.update(deltaTf(1:Nm));
                % 3. Update Rock Temperature
                Trm = ProductionSystem.Reservoir.State.Properties('Tr');
                Trm.update(deltaTr(1:Nm));
                
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
                        Tff = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('Tf');
                        Tff.update(deltaTf(IP:EP));
                        
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
        function [Tph, Gph, Tk] = TransmissibilityMatrix(obj, Grid, UpWind, Mob, rho, RhoInt)
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            % Transmissibility matrix construction
            Tx = zeros(Nx+1, Ny, Nz);
            Ty = zeros(Nx, Ny+1, Nz);
            Tz = zeros(Nx, Ny, Nz+1);
            THx = zeros(Nx+1, Ny, Nz);
            THy = zeros(Nx, Ny+1, Nz);
            THz = zeros(Nx, Ny, Nz+1);
            
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
            
            % Gravity Matrix
            Tx(2:Grid.Nx,:,:)= Tx(2:Grid.Nx,:,:) .* RhoInt.x(2:Grid.Nx,:,:);
            Ty(:,2:Grid.Ny,:)= Ty(:,2:Grid.Ny,:) .* RhoInt.y(:,2:Grid.Ny,:);
            Tz(:,:,2:Grid.Nz)= Tz(:,:,2:Grid.Nz) .* RhoInt.z(:,:,2:Grid.Nz);   
            Gph = ReshapeTransmisibility(Grid, Nx, Ny, Nz, N, Tx, Ty, Tz);
            
            % Rock Conductivity Matrix
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
            % 
        end
        function [Qf, Qhf, QTfTr]= ComputeSourceTerms_frac_mat(obj, ProductionSystem, DiscretizationModel)
            Qf = zeros(DiscretizationModel.N, obj.NofPhases);     % Mass Flow flux between each two media 
            Qhf = zeros(DiscretizationModel.N, obj.NofPhases);    % Heat Flow flux betweem each two media
            QTfTr = zeros(DiscretizationModel.N, obj.NofPhases);  % Convect-conduct heat flux between matrix rock and fracture fluid
            Nm = DiscretizationModel.ReservoirGrid.N;
            k_cond = ProductionSystem.Reservoir.k_cond(:,1);
%             U = zeros(DiscretizationModel.N,1);                   % Heat exchange coefficient

            % Global variables
            P = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'P_'); % useful for cross connections assembly
            rho = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'rho_'); % useful for cross connections assembly
            h = ProductionSystem.CreateGlobalVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], obj.NofPhases, 'h_'); % useful for cross connections assembly
            Tf = ProductionSystem.CreateGlobalSinglePhaseVariables([DiscretizationModel.ReservoirGrid, DiscretizationModel.FracturesGrid.Grids], 'Tf'); % useful for cross connections assembly
            Tr = ProductionSystem.Reservoir.State.Properties('Tr').Value;
            for ph=1:obj.NofPhases
                Index.Start = 1;
                Index.End = Nm;
                for c=1:length(DiscretizationModel.CrossConnections)
                    j = DiscretizationModel.CrossConnections(c).Cells;
                    i = c + Nm;
                    T_Geo = DiscretizationModel.CrossConnections(c).T_Geo;
                    UpWind = DiscretizationModel.CrossConnections(c).UpWind;
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
                     
                    QTfTr(i, ph)    = QTfTr(i, ph) + sum( T_Geo(j<=Nm) .* k_cond(j<=Nm) .* (Tr(j<=Nm, ph) - Tf(i, ph)) );
                    QTfTr(j<=Nm, ph) = QTfTr(j<=Nm, ph) + T_Geo(j<=Nm) .* k_cond(j<=Nm) .* (Tf(i, ph) - Tr(j<=Nm, ph));
                end
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
                [dQdp, dQdT] = Prod(i).dQdPdT(State, K, obj.Mob, obj.dMob, obj.drho, obj.NofPhases);
                [dQhdp, dQhdT] = Prod(i).dQhdPdT(State, K, obj.Mob, obj.dMob, obj.drho, obj.dh, obj.NofPhases);
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
        function ComputeTotalFluxes(obj, ProductionSystem, DiscretizationModel)
            % this is virtual call
        end
        function AverageMassOnCoarseBlocks(obj, ProductionSystem, FineGrid, FluidModel, ADMRest)
            T = ProductionSystem.CreateGlobalSinglePhaseVariables(FineGrid, 'Tf'); % No phase index for this variable
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
            T = ProductionSystem.Reservoir.State.Properties('Tf');
            rho = ProductionSystem.Reservoir.State.Properties('rho_1');
            T.update(delta_T(Start:End));
            rho.update(delta_rho(Start:End));
            
            % Updating the variables in fractures
            for frac = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                Start = End + 1;
                End = Start + FineGrid(frac+1).N - 1;
                T = ProductionSystem.FracturesNetwork.Fractures(frac).State.Properties('Tf');
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