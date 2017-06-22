% Immiscible Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 7 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Immiscible_formulation < formulation
    properties
        Mobt
    end
    methods
        function obj = Immiscible_formulation()
            obj@formulation();
            obj.Tph = cell(2,1);
            obj.Gph = cell(2,1);
        end
        function x = GetPrimaryUnknowns(obj, ProductionSystem, N)
             x = zeros(obj.NofPhases * N, 1);
             x(1:N) = ProductionSystem.Reservoir.State.Properties('P_2').Value;
             for i=1:obj.NofPhases-1
                 x(i*N + 1:(i+1)*N) = ProductionSystem.Reservoir.State.Properties(['S_', num2str(i)]).Value;
             end
        end
        function x = GetPrimaryPressure(obj, ProductionSystem, N)
            x = zeros(N, 1);
            if obj.NofPhases > 1
                x(1:N) = ProductionSystem.Reservoir.State.Properties('P_2').Value;
            else
                x(1:N) = ProductionSystem.Reservoir.State.Properties('P_1').Value;
            end
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            %% 1. Reservoir Properteis and Derivatives
            obj.drhodp = FluidModel.DrhoDp(ProductionSystem.Reservoir.State);
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.dPc = FluidModel.DPcDS(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            %% 2. Fractures Properteis and Derivatives
        end
        %% Methods for FIM Coupling
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Create local variables
            N = DiscretizationModel.ReservoirGrid.N;
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            s_old = zeros(N, obj.NofPhases);
            rho_old = zeros(N, obj.NofPhases);
            P = zeros(N, obj.NofPhases);
            s = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            
            % Copy values in local variables
            for i=1:obj.NofPhases
                s_old(:,i) = State0.Properties(['S_', num2str(i)]).Value;
                rho_old(:,i) = State0.Properties(['rho_', num2str(i)]).Value;
                P(:, i) = ProductionSystem.Reservoir.State.Properties(['P_', num2str(i)]).Value;
                s(:, i) = ProductionSystem.Reservoir.State.Properties(['S_', num2str(i)]).Value;
                rho(:, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            end
            depth = DiscretizationModel.ReservoirGrid.Depth;
            
            % Accumulation Term
            AS = speye(N)*pv/dt;
            
            % Transmissibility matrix
            for i=1:obj.NofPhases
                [obj.Tph{i}, obj.Gph{i}] = obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, obj.UpWind(i), obj.Mob(:,i), rho(:,i), obj.GravityModel.RhoInt(i));
            end
            % Source terms
            q = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            
            % RESIDUAL
            Residual = zeros(obj.NofPhases*N,1);
            for i=1:obj.NofPhases
                Residual((i-1)*N+1:i*N)  = AS*(rho(:,i) .* s(:,i) - rho_old(:,i) .* s_old(:,i))...
                                           + obj.Tph{i} * P(:, i)...
                                           - obj.Gph{i} * depth...
                                           - q(:,i);
            end   
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            % Create local variables
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N = DiscretizationModel.ReservoirGrid.N;
            pv = DiscretizationModel.ReservoirGrid.Volume*ProductionSystem.Reservoir.Por;
            rho = zeros(N, obj.NofPhases);
            s = zeros(N, obj.NofPhases);
            for i=1:obj.NofPhases
                rho(:,i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
                s(:,i) = ProductionSystem.Reservoir.State.Properties(['S_', num2str(i)]).Value;
            end
            
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            Jp = cell(2,1);
            JS = cell(2,1);
            
            for i=1:2
                % 1.a Pressure Block
                Jp{i} = obj.Tph{i};

                % 1.b: compressibility part
                dMupx = obj.UpWind(i).x*(obj.Mob(:, i) .* obj.drhodp(:,i));
                dMupy = obj.UpWind(i).y*(obj.Mob(:, i) .* obj.drhodp(:,i));
                dMupz = obj.UpWind(i).z*(obj.Mob(:, i) .* obj.drhodp(:,i));
                
                vecX1 = min(reshape(obj.U(i).x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                vecX2 = max(reshape(obj.U(i).x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                vecY1 = min(reshape(obj.U(i).y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                vecY2 = max(reshape(obj.U(i).y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                vecZ1 = min(reshape(obj.U(i).z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                vecZ2 = max(reshape(obj.U(i).z(:,:,2:Nz+1), N, 1), 0) .* dMupz; 
                acc = pv/dt .* obj.drhodp(:,i) .* s(:,i);
                
                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                Jp{i} = Jp{i} + spdiags(DiagVecs, DiagIndx, N, N);
                
                % 2. Saturation Block
                dMupx = obj.UpWind(i).x * (obj.dMob(:,i) .* rho(:,i));
                dMupy = obj.UpWind(i).y * (obj.dMob(:,i) .* rho(:,i));
                dMupz = obj.UpWind(i).z * (obj.dMob(:,i) .* rho(:,i));
                % Construct JS block
                x1 = min(reshape(obj.U(i).x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                x2 = max(reshape(obj.U(i).x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                y1 = min(reshape(obj.U(i).y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                y2 = max(reshape(obj.U(i).y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                z1 = min(reshape(obj.U(i).z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                z2 = max(reshape(obj.U(i).z(:,:,2:Nz+1), N, 1), 0) .* dMupz;
                v = (-1)^(i+1) * ones(N,1)*pv/dt .* rho(:,i);
                DiagVecs = [-z2, -y2, -x2, z2+y2+x2-z1-y1-x1+v, x1, y1, z1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                JS{i} = spdiags(DiagVecs,DiagIndx,N,N);
            end  
            
            %Add capillarity
            JS {1}= JS{1} - Jp{1} * spdiags(obj.dPc, 0, N, N);
            
            % Add Wells
            [Jp{1}, JS{1}, Jp{2}, JS{2}] = obj.AddWellsToJacobian(Jp{1}, JS{1}, Jp{2}, JS{2}, ProductionSystem.Reservoir.State, ProductionSystem.Wells, ProductionSystem.Reservoir.K(:,1));
            
            % Full Jacobian: put the 4 blocks together
            Jacobian = [Jp{1}, JS{1}; Jp{2}, JS{2}];
        end
        function delta = UpdateState(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                return
            else
                Nm =  DiscretizationModel.ReservoirGrid.N;
                %% 1. Update matrix
                % Update Pressure
                Pm = ProductionSystem.Reservoir.State.Properties(['P_', num2str(obj.NofPhases)]);
                Pm.update(delta(1:Nm));
                DeltaLast = zeros(Nm, 1);
                for ph = 1:obj.NofPhases-1
                    Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(ph)]);
                    DeltaS = delta(ph*Nm + 1:(ph+1)*Nm);
                    Sm.update(DeltaS);
                    % Remove values that are not physical
                    Sm.Value = max(Sm.Value, 0);
                    Sm.Value = min(Sm.Value, 1);
                    DeltaLast = DeltaLast + DeltaS;
                end
                Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(obj.NofPhases)]);
                Sm.update(-DeltaLast);
                % Remove values that are not physical
                Sm.Value = max(Sm.Value, 0);
                Sm.Value = min(Sm.Value, 1);
                % Update Phase Densities
                FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
                % Update total density
                FluidModel.ComputeTotalDensity(ProductionSystem.Reservoir.State);
                % Update Pc
                FluidModel.ComputePc(ProductionSystem.Reservoir.State);
                
                %% 2. Update fractures pressure and densities
                if ProductionSystem.FracturesNetwork.Active
                    for i=1:ProductionSystem.FracturesNetwork.NofFractures
                        % Update Pressure
                        Pf = ProductionSystem.FracturesNetwork.Fractures(i).State.Properties(['P_', num2str(obj.NofPhases)]);
                        Pf.update(delta);
                        DeltaLast = zeros(Nf(i), 1);
                        for ph = 1:obj.NofPhases-1
                            Sf = ProductionSystem.FracturesNetwork.Fractures(i).State.Properties(['S_', num2str(ph)]);
                            Sf.update(delta(ph*Nf(i) + 1:(ph+1)*Nf(i)));
                            Sf.Value = max(Sf.Value, 0);
                            Sf.Value = min(Sf.Value, 1);
                            DeltaLast = DeltaLast + delta(ph*Nf(i) + 1:(ph+1)*Nf(i));
                        end
                        Sf = ProductionSystem.FracturesNetwork.Fractures(i).State.Properties(['S_', num2str(obj.NofPhases)]);
                        Sf.update(-DeltaLast);
                        Sf.Value = max(Sf.Value, 0);
                        Sf.Value = min(Sf.Value, 1);
                        % Update Phase Densities
                        FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(i).State);
                        % Update total density
                        FluidModel.ComputeTotalDensity(ProductionSystem.FracturesNetwork.Fractures(i).State);
                        % Update Pc
                        FluidModel.ComputePc(ProductionSystem.FracturesNetwork.Fractures(i).State);
                    end
                end
            end
        end
        function [Tph, Gph] = TransmissibilityMatrix(obj, Grid, UpWind, Mob, rho, RhoInt)
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            % Transmissibility matrix construction
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
            % Construct matrix
            x1 = reshape(Tx(1:Nx,:,:), N, 1);
            x2 = reshape(Tx(2:Nx+1,:,:), N, 1);
            y1 = reshape(Ty(:,1:Ny,:), N, 1);
            y2 = reshape(Ty(:,2:Ny+1,:), N, 1);
            z1 = reshape(Tz(:,:,1:Nz), N, 1);
            z2 = reshape(Tz(:,:,2:Nz+1), N, 1);
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
            Tph = spdiags(DiagVecs,DiagIndx,N,N);
            
            % Gravity Matrix
            Tx(2:Grid.Nx,:,:)= Tx(2:Grid.Nx,:,:) .* RhoInt.x(2:Grid.Nx,:,:);
            Ty(:,2:Grid.Ny,:)= Ty(:,2:Grid.Ny,:) .* RhoInt.y(:,2:Grid.Ny,:);
            Tz(:,:,2:Grid.Nz)= Tz(:,:,2:Grid.Nz) .* RhoInt.z(:,:,2:Grid.Nz);
            
            % Construct matrix
            x1 = reshape(Tx(1:Nx,:,:), N, 1);
            x2 = reshape(Tx(2:Nx+1,:,:), N, 1);
            y1 = reshape(Ty(:,1:Ny,:), N, 1);
            y2 = reshape(Ty(:,2:Ny+1,:), N, 1);
            z1 = reshape(Tz(:,:,1:Nz), N, 1);
            z2 = reshape(Tz(:,:,2:Nz+1), N, 1);
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
            Gph = spdiags(DiagVecs, DiagIndx, N, N);
        end
        function q = ComputeSourceTerms(obj, N, Wells)
            q = zeros(N, obj.NofPhases);    
            %Injectors
            for i=1:Wells.NofInj
                c = Wells.Inj(i).Cells;
                q(c, :) = Wells.Inj(i).QPhases(:,:);
            end
            %Producers
            for i=1:Wells.NofProd
                c = Wells.Prod(i).Cells;
                q(c, :) = Wells.Prod(i).QPhases(:,:);
            end
        end
        function [Jwp, JwS, Jnwp, JnwS] = AddWellsToJacobian(obj, Jwp, JwS, Jnwp, JnwS, State, Wells, K)
            % Define Local handles
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            p = State.Properties('P_2').Value;
            rho = zeros(length(p), obj.NofPhases);
            for i = 1:obj.NofPhases
                rho(:, i) = State.Properties(['rho_', num2str(i)]).Value;
            end
            
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                for j=1:length(a)
                    Jnwp(a(j),a(j)) = Jnwp(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mob(:, 2)*Inj(i).rho(j, 2);
                    Jwp(a(j),a(j)) = Jwp(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mob(:, 1)*Inj(i).rho(j, 1);
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                for j=1:length(b)
                    Jnwp(b(j),b(j)) = Jnwp(b(j),b(j)) + Prod(i).PI*K(b(j)).*obj.Mob(b(j), 2) .* rho(b(j), 2)...
                     - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * obj.drhodp(b(j),2) .* (Prod(i).p(j) - p(b(j)));                    
                    Jwp(b(j),b(j)) = Jwp(b(j),b(j)) + Prod(i).PI*K(b(j)).*obj.Mob(b(j), 1) .* rho(b(j), 1)...
                     - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) .* obj.drhodp(b(j),1) .* (Prod(i).p(j) - p(b(j)));
                    
                    JwS(b(j),b(j)) = JwS(b(j),b(j)) - Prod(i).PI*K(b(j)).* rho(b(j), 1) .* (Prod(i).p(j) - p(b(j))).*obj.dMob(b(j), 1);
                    JnwS(b(j),b(j)) = JnwS(b(j),b(j)) - Prod(i).PI*K(b(j)).* rho(b(j), 2) .* (Prod(i).p(j) - p(b(j))).*obj.dMob(b(j), 2);
                end
            end
        end
        function AverageMassOnCoarseBlocks(obj, Status, FluidModel, R, P)
            % Perform Average for ADM
            S_rest = R * Status.S;
            Status.S = P * S_rest;
            % Update other unknwons as well 
            %obj.UpdatePhaseCompositions(Status, FluidModel);
        end
        function CFL = ComputeCFLNumber(obj, ProductionSystem, DiscretizationModel, dt)
            N = DiscretizationModel.ReservoirGrid.N;      
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            P = zeros(N, obj.NofPhases);
            S = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            
            % Copy values in local variables
            for i=1:obj.NofPhases
                P(:, i) = ProductionSystem.Reservoir.State.Properties(['P_', num2str(i)]).Value;
                rho(:, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
                S(:, i) = ProductionSystem.Reservoir.State.Properties(['S_', num2str(i)]).Value;
            end
            
            for i=1:obj.NofPhases
                [obj.Tph{i}, obj.Gph{i}] = obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, obj.UpWind(i), obj.Mob(:,i), rho(:,i), obj.GravityModel.RhoInt(i));
            end
             % Depths
            depth = DiscretizationModel.ReservoirGrid.Depth;
            
            % Source terms
            q = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            
            ThroughPut = zeros(N, obj.NofPhases);
            Mass = zeros(N, obj.NofPhases);
            for i=1:obj.NofPhases
                % extract lower diagonals
                d = tril(obj.Tph{i}, -1);
                % extract upper diagonals
                u = triu(obj.Tph{i},  1);
                % extract main diagonal
                Diag = diag(obj.Tph{1});
                % Assemble
                D = diag(Diag + sum(u, 2)) + d; 
                U = diag(Diag + sum(d, 2)) + u;
                % extract lower diagonals
                gd = tril(obj.Gph{i}, -1);
                % extract upper diagonals
                gu = triu(obj.Gph{i},  1);
                % extract main diagonal
                DiagG = diag(obj.Gph{1});
                % Assemble
                GD = diag(DiagG + sum(gu, 2)) + gd; 
                GU = diag(DiagG + sum(gd, 2)) + gu;
                ThroughPut(:,i) = ...
                           - min(D * P(:,i), 0)...    % Convective term (take only incoming fluxes (negative))                
                           - min(U * P(:,i), 0)...    % Convective term (take only incoming fluxes (negative))     
                           - min(GD * depth, 0)...    % Gravity term (take only incoming fluxes (negative))                
                           - min(GU * depth, 0)...    % Gravity term (take only incoming fluxes (negative))
                           + max(q(:,i), 0);          % Wells (injectors)
                Mass(:,i) = rho(:,i) .* S(:,i);
            end
            Mass = max(Mass, 1e-10);
            ThroughPut(ThroughPut < 1e-10) = 0;
            Ratio = ThroughPut ./ Mass;
            CFL = dt/pv * max(max(Ratio));
        end
        %% Methods for Sequential Coupling
        function ComputeTotalMobility(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.Mobt = sum(obj.Mob, 2);
        end
        function UpdateFractionalFlow(obj, ProductionSystem, FluidModel)
            obj.ComputeTotalMobility(ProductionSystem, FluidModel);
            obj.f = obj.Mob(:,1) ./ obj.Mobt;
        end
        function dfdS(obj, ProductionSystem, FluidModel)
            dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State.S);
            obj.df = (dMob(:,1) .* sum(obj.Mob, 2) - sum(dMob, 2) .* obj.Mob(:,1)) ./ sum(obj.Mob, 2).^2;
        end
        function Residual = BuildPressureResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Create local variables
            N = DiscretizationModel.ReservoirGrid.N;
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            s_old = zeros(N, obj.NofPhases);
            rho_old = zeros(N, obj.NofPhases);
            P = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            
            % Copy values in local variables
            for i=1:obj.NofPhases
                P(:, i) = ProductionSystem.Reservoir.State.Properties(['P_', num2str(i)]).Value;
                s_old(:, i) = State0.Properties(['S_', num2str(i)]).Value;
                rho_old(:, i) = State0.Properties(['rho_', num2str(i)]).Value;
                rho(:, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            end 
            depth = DiscretizationModel.ReservoirGrid.Depth;
            
            
            % Transmissibility matrix
            for i=1:obj.NofPhases
                [obj.Tph{i}, ~] = obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, obj.UpWind(i), obj.Mob(:,i), rho(:,i), obj.GravityModel.RhoInt(i));
                %obj.Tph{i} = bsxfun(@rdivide, obj.Tph{i}, rho(:, i));
            end
            
            % Source terms
            q = sparse(obj.ComputeSourceTerms(N, ProductionSystem.Wells));
            
            %% RESIDUAL
            Residual = zeros(N, 1);
            for i=1:obj.NofPhases
            Residual(:) = ...
                    Residual(:) - ...
                    pv/dt * (rho_old(:, i) .* s_old(:, i) ./ rho(:, i))...
                    + obj.Tph{i} * P(:, i)...
                    - q(:, i) ./ rho(:, i);
            end
            Residual(:) = Residual(:) + pv/dt;
            
            Residual = zeros(N, 1);
            for i=1:obj.NofPhases
            Residual(:) = ...
                   Residual(:) + ...
                   pv/dt * (rho(:, i).*s_old(:, i) - rho_old(:, i) .* s_old(:, i))...
                   + obj.Tph{i} * P(:, i)...
                   - q(:, i);
            end
        end
        function J = BuildPressureMatrix(obj, ProductionSystem, DiscretizationModel, dt, State0)
            %% 0. Initialise local variables
            N = DiscretizationModel.ReservoirGrid.N;
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            s_old = zeros(N, obj.NofPhases);
            rho_old = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            
            % Copy values in local variables
            for i=1:obj.NofPhases
                s_old(:, i) = State0.Properties(['S_', num2str(i)]).Value;
                rho_old(:, i) = State0.Properties(['rho_', num2str(i)]).Value;
                rho(:, i) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value;
            end 
            
            A = sparse(N,N);
            C = sparse(N,N);
            for i=1:obj.NofPhases
                %% 1. Accummulation term
                vec = - pv/dt * (- rho_old(:, i) .* s_old(:, i) .* obj.drhodp(:,i) ./ rho(:, i).^2) ;
                % Alternative
                vec = pv/dt * obj.drhodp(:,i);
                C = C + spdiags(vec, 0, N, N);
                
                %% 2. Convective term 
                % a. transmissibility matrix
                A = A + obj.Tph{i};
                % b. Derivative of transmissibility with respect to p
                dMupx = obj.UpWind(i).x*(obj.Mob(:, i) .* obj.drhodp(:,i));
                dMupy = obj.UpWind(i).y*(obj.Mob(:, i) .* obj.drhodp(:,i));
                dMupz = obj.UpWind(i).z*(obj.Mob(:, i) .* obj.drhodp(:,i));
                
                vecX1 = min(reshape(obj.U(i).x(1:Nx,:,:), N, 1), 0)   .* dMupx;
                vecX2 = max(reshape(obj.U(i).x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
                vecY1 = min(reshape(obj.U(i).y(:,1:Ny,:), N, 1), 0)   .* dMupy;
                vecY2 = max(reshape(obj.U(i).y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
                vecZ1 = min(reshape(obj.U(i).z(:,:,1:Nz), N, 1), 0)   .* dMupz;
                vecZ2 = max(reshape(obj.U(i).z(:,:,2:Nz+1), N, 1), 0) .* dMupz; 
                maindiag = vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1;
                
                DiagVecs = [-vecZ2, -vecY2, -vecX2, maindiag, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                dTdp = spdiags(DiagVecs, DiagIndx, N, N);
%                 
%                 % multiply row i by rho(i)
%                 dTdp = bsxfun(@times, dTdp, rho(:, i));
%                 
%                 dMupx = obj.UpWind(i).x*(obj.Mob(:, i) .* rho(:,i));
%                 dMupy = obj.UpWind(i).y*(obj.Mob(:, i) .* rho(:,i));
%                 dMupz = obj.UpWind(i).z*(obj.Mob(:, i) .* rho(:,i));
%                 
%                 vecX1 = min(reshape(obj.U(i).x(1:Nx,:,:), N, 1), 0)   .* dMupx;
%                 vecX2 = max(reshape(obj.U(i).x(2:Nx+1,:,:), N, 1), 0) .* dMupx;
%                 vecY1 = min(reshape(obj.U(i).y(:,1:Ny,:), N, 1), 0)   .* dMupy;
%                 vecY2 = max(reshape(obj.U(i).y(:,2:Ny+1,:), N, 1), 0) .* dMupy;
%                 vecZ1 = min(reshape(obj.U(i).z(:,:,1:Nz), N, 1), 0)   .* dMupz;
%                 vecZ2 = max(reshape(obj.U(i).z(:,:,2:Nz+1), N, 1), 0) .* dMupz; 
%                 maindiag = vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1;
%                 vec = - maindiag .* obj.drhodp(:,i); 
%                 dTdp = dTdp + spdiags(vec, 0, N, N);
%                  
%                 dTdp = bsxfun(@rdivide, dTdp, rho(:, i).^2);
%                  
                 A = A + dTdp; 
            end
            
            %% 3. Add wells
            W = obj.AddWellsToPressureSystem(N, ProductionSystem.Reservoir.State, ProductionSystem.Wells, ProductionSystem.Reservoir.K(:,1), rho);
            
            %% Jacobian: Put them together
            J = A + C + W;
        end       
        function W = AddWellsToPressureSystem(obj, N, State, Wells, K, rho)
            %% Add Wells in residual form
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            dq = zeros(N, obj.NofPhases);
            %Injectors
            for i=1:length(Inj)
                c = Inj(i).Cells;
                [dq(c, :), ~] = Inj(i).dQPhasesdPdS(State, K, obj.NofPhases);
            end
            
            %Producers
            for i=1:length(Prod)
                c = Prod(i).Cells;
                [dq(c, :), ~] = Prod(i).dQPhasesdPdS(State, K, obj.Mob, obj.dMob, obj.drhodp, obj.NofPhases);
            end
            q = sparse(obj.ComputeSourceTerms(N, Wells));
            dqt = sum((dq .* rho - obj.drhodp .* q)./rho.^2, 2);
            % Alternative
            dqt = sum(dq, 2);
            
            W = -spdiags(dqt, 0, N, N);
            
        end
        function UpdatePressure(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            %% 1. Update matrix pressure and densities
            % Update Pressure
            Pm = ProductionSystem.Reservoir.State.Properties(['P_', num2str(obj.NofPhases)]);
            Pm.update(delta);
            % Update Phase Densities
            FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
            % Update total density
            FluidModel.ComputeTotalDensity(ProductionSystem.Reservoir.State);
            if obj.NofPhases > 1
            % Update Pc
                FluidModel.ComputePc(ProductionSystem.Reservoir.State);
            end
            
            %% 2. Update fractures pressure and densities
            if ProductionSystem.FracturesNetwork.Active
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    % Update Pressure
                    Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_', num2str(obj.NofPhases)]);
                    Pf.update(delta);
                    % Update Phase Densities
                    FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    % Update total density
                    FluidModel.ComputeTotalDensity(ProductionSystem.FracturesNetwork.Fractures(f).State);
                end
            end
        end
        function [Utot, Qwells] = ComputeTotalFluxes(obj, ProductionSystem, DiscretizationModel)
            Utot.x(2:Nx+1,:,:) = obj.U(2).x(2:Nx+1,:,:) .* reshape(obj.UpWind.x *  obj.Mobt, Nx, Ny, Nz); %- Ucap.x(2:Nx,:);
            Utot.y(:,2:Ny+1,:) = obj.U(2).y(:,2:Ny+1,:) .* reshape(obj.UpWind.y *  obj.Mobt, Nx, Ny, Nz); %- Ucap.y(:,2:Ny);
            Utot.z(:,:,2:Nz+1) = obj.U(2).z(:,:,2:Nz+1) .* reshape(obj.UpWind.z *  obj.Mobt, Nx, Ny, Nz);  %- Ucap.y(:,2:Ny);
            
            % Wells total fluxes
            Qwells = ProductionSystem.Wells.TotalFluxes(ProductionSystem.Reservoir, obj.Mobt);
        end
        function conservative = CheckMassConservation(obj, Grid, Utot, Qwells)
            %Checks mass balance in all cells
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            conservative = 1;
            maxUx = max(max(max(obj.Utot.x)));
            maxUy = max(max(max(obj.Utot.y)));
            maxUz = max(max(max(obj.Utot.z)));
            maxU = max([maxUx, maxUy, maxUz]);
            qWells = reshape(Qwells, Nx, Ny, Nz);
            for k=1:Nz
                for j=1:Ny
                    for i=1:Nx
                        Accum = Utot.x(i,j,k) - Utot.x(i+1,j,k) + Utot.y(i,j,k) - Utot.y(i,j+1,k) + Utot.z(i,j,k) - Utot.z(i,j,k+1) + qWells(i,j,k);
                        if (abs(Accum/maxU) > 10^(-5))
                            conservative = 0;
                        end
                    end
                end
            end
        end
        function ViscousMatrix(obj, Grid, Utot)
            %Builds Upwind Flux matrix
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;                                   
            q = min(obj.Qwells, 0);                        
            % right to left and top to bottom (negative x, y, z)
            Xneg = min(Utot.x, 0); 
            Yneg = min(Utot.y, 0);
            Zneg = min(Utot.z, 0);
            % make them vectors 
            x1 = reshape(Xneg(1:Nx,:,:),N,1);
            y1 = reshape(Yneg(:,1:Ny,:),N,1);
            z1 = reshape(Zneg(:,:,1:Nz),N,1);
            
            % left to right and bottom to top (positive x, y, z)
            Xpos = max(Utot.x, 0); 
            Ypos = max(Utot.y, 0); 
            Zpos = max(Utot.z, 0);
            % make them vectors
            x2 = reshape(Xpos(2:Nx+1,:,:), N, 1);
            y2 = reshape(Ypos(:,2:Ny+1,:), N, 1);
            z2 = reshape(Zpos(:,:,2:Nz+1), N, 1);
            
            % Assemble matrix
            DiagVecs = [z2, y2, x2, q+x1-x2+y1-y2+z1-z2, -x1, -y1, -z1]; % diagonal vectors
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny]; % diagonal index
            obj.V = spdiags(DiagVecs, DiagIndx, N, N);
        end
        function Residual = BuildTransportResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Initialise local objects
            pv = ProductionSystem.Reservoir.Por * DiscretizationModel.ReservoirGrid.Volume;
            s = ProductionSystem.Reservoir.State.S;
            s_old = State0.S;      
            
            % Compute residual
            Residual = pv/dt * (s - s_old)  - max(obj.Qwells, 0) - obj.V * obj.f;
        end
        function Jacobian = BuildTransportJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            % Build Transport Jacobian
            pv = ProductionSystem.Reservoir.Por * DiscretizationModel.ReservoirGrid.Volume;
            N = DiscretizationModel.ReservoirGrid.N;
            D = spdiags(pv/dt*ones(N,1),0,N,N);
            Jacobian = D - obj.V * spdiags(obj.df,0,N,N); %+ CapJac;
        end
        function UpdateSaturation(obj, State, delta, FluidModel)
            State.S = State.S + delta;
            State.S = min(State.S, 1);
            State.S = max(State.S, FluidModel.Phases(1).sr);
        end
        function UpdateSaturationExplicitly(obj, ProductionSystem, DiscretizationModel, dt)
            % 0. Initialise
            pv = ProductionSystem.Reservoir.Por .* DiscretizationModel.ReservoirGrid.Volume;
            N = DiscretizationModel.ReservoirGrid.N;
            
            % 1. Solve
            T = spdiags(dt/pv*ones(N,1),0,N,N);    % dt/pv * Cell Fluxes and producer
            B = T * obj.V;
            injector = max(obj.Qwells,0) .* dt/pv;  % injection flux * dt/pv
            
            ProductionSystem.Reservoir.State.S = ProductionSystem.Reservoir.State.S + (B * obj.f + injector);
        end
    end
end