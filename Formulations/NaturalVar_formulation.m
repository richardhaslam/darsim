% Natural variable Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef NaturalVar_formulation < Compositional_formulation
    properties
        K
        dKdp
        InitialPhaseState
        PreviousSinglePhase
    end
    methods
        function obj = NaturalVar_formulation(n_cells, n_components)
            obj@Compositional_formulation(n_components);
            obj.PreviousSinglePhase = zeros(n_cells, 1);
            obj.SinglePhase = zeros(n_cells, 1);
        end
        function SavePhaseState(obj)
            obj.InitialPhaseState = obj.SinglePhase;
        end
        function Reset(obj)
            obj.SinglePhase = obj.InitialPhaseState;
            obj.PreviousSinglePhase = obj.InitialPhaseState;
        end
        function x = GetPrimaryUnknowns(obj, ProductionSystem, DiscretizationModel)
            N =DiscretizationModel.N;
            x = zeros(obj.NofComponents * N, 1);
            x(1:N) = ProductionSystem.Reservoir.State.Properties('P_2').Value;
            for i=1:obj.NofPhases
                x(i*N + 1:(i+1)*N) = ProductionSystem.Reservoir.State.Properties(['S_', num2str(i)]).Value;
            end
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.drhodp = FluidModel.DrhoDp(ProductionSystem.Reservoir.State, obj.SinglePhase);
            obj.dPc = FluidModel.DPcDS(ProductionSystem.Reservoir.State.Properties('S_1').Value);
            obj.K = FluidModel.ComputeKvalues(ProductionSystem.Reservoir.State);
            obj.dKdp = FluidModel.DKvalDp(ProductionSystem.Reservoir.State);
            obj.SinglePhase = FluidModel.CheckNumberOfPhases(obj.SinglePhase, obj.PreviousSinglePhase, ProductionSystem.Reservoir.State, obj.K);
         end
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)                   
            %Create local variables
            N = DiscretizationModel.ReservoirGrid.N;          
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            S_old = zeros(N, obj.NofPhases);
            rho_old = zeros(N, obj.NofPhases);
            x_old = zeros(N, obj.NofComponents*obj.NofPhases);
            P = zeros(N, obj.NofPhases);
            S = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            x = zeros(N, obj.NofComponents*obj.NofPhases);
            
            % Copy values in local variables
            for j=1:obj.NofPhases
                rho_old(:,j) = State0.Properties(['rho_', num2str(j)]).Value;
                P(:, j) = ProductionSystem.Reservoir.State.Properties(['P_', num2str(j)]).Value;
                rho(:, j) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(j)]).Value;
                S_old(:,j) = State0.Properties(['S_', num2str(j)]).Value;
                S(:, j) = ProductionSystem.Reservoir.State.Properties(['S_', num2str(j)]).Value;
                for i=1:obj.NofComponents
                    x_old(:,(i-1)*obj.NofPhases + j) = State0.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
                    x(:,(i-1)*obj.NofPhases + j) = ProductionSystem.Reservoir.State.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
                end
            end
            
            % Depths
            depth = DiscretizationModel.ReservoirGrid.Depth;
            
            %Accumulation term
            A = speye(N)*pv/dt;
                     
            %Source terms
            q = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            
            %% RESIDUAL
            Rbalance = zeros(N*obj.NofComponents, 1);
            Req = zeros(N*obj.NofComponents, 1);
            for i=1:obj.NofComponents
                % 1. MASS CONSERVATION EQUATIONS                 
                m = x(:,(i-1)*2+1) .* rho(:,1) .* S(:,1) + x(:,(i-1)*2+2) .* rho(:,2) .* S(:,2);
                m_old = x_old(:,(i-1)*2+1) .* rho_old(:,1) .* S_old(:,1) + x_old(:,(i-1)*2+2) .* rho_old(:,2) .* S_old(:,2);
                % Phase Transmissibilities
                obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, rho, obj.GravityModel.RhoInt, x(:,(i-1)*2+1:(i-1)*2+2), i);          
                % Residual
                Rbalance((i-1)*N+1:i*N) = ...
                    A * m - A * m_old...          % Accumulation term
                    + obj.Tph{i, 1} *  P(:,1) ... % Convective term                
                    + obj.Tph{i, 2} *  P(:,2)...
                    - obj.Gph{i,1} * depth...     % Gravity
                    - obj.Gph{i,2} * depth...
                    - q(:,i);                     % Source terms
                
                % 2. THERMODYNAMIC EQUILIBRIUM EQUATIONS
                Rcomp = x(:,(i-1)*2 + 1) - obj.K(:,i).*x(:,(i-1)*2 + 2);
                % If Single phase set residual to be zero
                Rcomp(obj.SinglePhase > 0) = 0;
                Req((i-1)*N+1:i*N) = Rcomp;                 
            end
            
            %Stick them together
            Residual = [Rbalance; Req];
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            
            % Initialise local variables
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N = DiscretizationModel.ReservoirGrid.N;
            pv = DiscretizationModel.ReservoirGrid.Volume*ProductionSystem.Reservoir.Por;
            
            P = zeros(N, obj.NofPhases);
            rho = zeros(N, obj.NofPhases);
            S = zeros(N, obj.NofPhases);
            x = zeros(N, obj.NofComponents*obj.NofPhases);
            
            % Copy values in local variables
            for j=1:obj.NofPhases
                P(:, j) = ProductionSystem.Reservoir.State.Properties(['P_', num2str(j)]).Value;
                rho(:, j) = ProductionSystem.Reservoir.State.Properties(['rho_', num2str(j)]).Value;
                S(:, j) = ProductionSystem.Reservoir.State.Properties(['S_', num2str(j)]).Value;
                for i=1:obj.NofComponents
                    x(:,(i-1)*obj.NofPhases + j) = ProductionSystem.Reservoir.State.Properties(['x_', num2str(i),'ph',num2str(j)]).Value;
                end
            end
            
            % Fill in block by block
            Jp = cell(obj.NofComponents, 1);
            JS = cell(obj.NofComponents, 1);
            for i=1:obj.NofComponents
                %% 1. Component i pressure block
                % 1.a: divergence
                Jp{i} = obj.Tph{i,1}  + obj.Tph{i, 2};
                % 1.b: compressibility part
                dMupxPh1 = obj.UpWind{1}.x * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drhodp(:,1));
                dMupyPh1 = obj.UpWind{1}.y * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drhodp(:,1));
                dMupzPh1 = obj.UpWind{1}.z * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drhodp(:,1));
                dMupxPh2 = obj.UpWind{2}.x * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drhodp(:,2));
                dMupyPh2 = obj.UpWind{2}.y * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drhodp(:,2));
                dMupzPh2 = obj.UpWind{2}.z * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drhodp(:,2));
                
                vecX1 = min(reshape(obj.U{1}.x(1:Nx,:,:),N,1), 0).*dMupxPh1 + min(reshape(obj.U{2}.x(1:Nx,:,:),N,1), 0).*dMupxPh2;
                vecX2 = max(reshape(obj.U{1}.x(2:Nx+1,:,:),N,1), 0).*dMupxPh1 + max(reshape(obj.U{2}.x(2:Nx+1,:,:),N,1), 0).*dMupxPh2;
                vecY1 = min(reshape(obj.U{1}.y(:,1:Ny,:),N,1), 0).*dMupyPh1 + min(reshape(obj.U{2}.y(:,1:Ny,:),N,1), 0).*dMupyPh2;
                vecY2 = max(reshape(obj.U{1}.y(:,2:Ny+1,:),N,1), 0).*dMupyPh1 + max(reshape(obj.U{2}.y(:,2:Ny+1,:),N,1), 0).*dMupyPh2;
                vecZ1 = min(reshape(obj.U{1}.z(:,:,1:Nz),N,1), 0).*dMupzPh1 + min(reshape(obj.U{2}.z(:,:,1:Nz),N,1), 0).*dMupzPh2;
                vecZ2 = max(reshape(obj.U{1}.z(:,:,2:Nz+1),N,1), 0).*dMupzPh1 + max(reshape(obj.U{2}.z(:,:,2:Nz+1),N,1), 0).*dMupzPh2;
                
                acc = pv/dt .* ( x(:,(i-1)*2+1) .* obj.drhodp(:,1) .* S(:,1) + x(:,(i-1)*2+2) .* obj.drhodp(:,2) .* S(:,2));
                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, -vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                Jp{i} = Jp{i} + spdiags(DiagVecs, DiagIndx, N, N);
                
                %% 2. Component i saturation block
                dMupxPh1 = obj.UpWind{1}.x*(obj.dMob(:,1) .* x(:,(i-1)*2+1) .* rho(:,1));
                dMupyPh1 = obj.UpWind{1}.y*(obj.dMob(:,1) .* x(:,(i-1)*2+1) .* rho(:,1));
                dMupzPh1 = obj.UpWind{1}.z*(obj.dMob(:,1) .* x(:,(i-1)*2+1) .* rho(:,1));
                dMupxPh2 = obj.UpWind{2}.x*(obj.dMob(:,2) .* x(:,(i-1)*2+2) .* rho(:,2));
                dMupyPh2 = obj.UpWind{2}.y*(obj.dMob(:,2) .* x(:,(i-1)*2+2) .* rho(:,2));
                dMupzPh2 = obj.UpWind{2}.z*(obj.dMob(:,2) .* x(:,(i-1)*2+2) .* rho(:,2));
                
                vecX1 = min(reshape(obj.U{1}.x(1:Nx,:,:),N,1), 0).*dMupxPh1 + min(reshape(obj.U{2}.x(1:Nx,:,:),N,1), 0).*dMupxPh2;
                vecX2 = max(reshape(obj.U{1}.x(2:Nx+1,:,:),N,1), 0).*dMupxPh1 + max(reshape(obj.U{2}.x(2:Nx+1,:,:),N,1), 0).*dMupxPh2;
                vecY1 = min(reshape(obj.U{1}.y(:,1:Ny,:),N,1), 0).*dMupyPh1 + min(reshape(obj.U{2}.y(:,1:Ny,:),N,1), 0).*dMupyPh2;
                vecY2 = max(reshape(obj.U{1}.y(:,2:Ny+1,:),N,1), 0).*dMupyPh1 + max(reshape(obj.U{2}.y(:,2:Ny+1,:),N,1), 0).*dMupyPh2;
                vecZ1 = min(reshape(obj.U{1}.z(:,:,1:Nz),N,1), 0).*dMupzPh1 + min(reshape(obj.U{2}.z(:,:,1:Nz),N,1), 0).*dMupzPh2;
                vecZ2 = max(reshape(obj.U{1}.z(:,:,2:Nz+1),N,1), 0).*dMupzPh1 + max(reshape(obj.U{2}.z(:,:,2:Nz+1),N,1), 0).*dMupzPh2;

                acc = pv/dt .* (x(:,(i-1)*2+1) .* rho(:,1) - x(:,(i-1)*2+2) .* rho(:,2));
                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                JS{i} = spdiags(DiagVecs,DiagIndx, N, N);
                % Capillarity
                JS{i} = JS{i} - obj.Tph{i,1} * spdiags(obj.dPc, 0, N, N);               
            end
            
            %% 3. Component 1 x1ph1 block
            dMupxPh1 = obj.UpWind{1}.x * (obj.Mob(:,1) .* rho(:,1));
            dMupyPh1 = obj.UpWind{1}.y * (obj.Mob(:,1) .* rho(:,1));
            dMupzPh1 = obj.UpWind{1}.z * (obj.Mob(:,1) .* rho(:,1));
            vecX1 = min(reshape(obj.U{1}.x(1:Nx,:,:),N,1), 0).*dMupxPh1; 
            vecX2 = max(reshape(obj.U{1}.x(2:Nx+1,:,:),N,1), 0).*dMupxPh1; 
            vecY1 = min(reshape(obj.U{1}.y(:,1:Ny,:),N,1), 0).*dMupyPh1; 
            vecY2 = max(reshape(obj.U{1}.y(:,2:Ny+1,:),N,1), 0).*dMupyPh1;
            vecZ1 = min(reshape(obj.U{1}.z(:,:,1:Nz),N,1), 0).*dMupzPh1; 
            vecZ2 = max(reshape(obj.U{1}.z(:,:,2:Nz+1),N,1), 0).*dMupzPh1;
            acc = pv/dt .* (S(:,1) .* rho(:,1));
            
            DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J1x1ph1 = spdiags(DiagVecs,DiagIndx, N, N);
            
            %% 4. Component 1 x1ph2  block
            dMupxPh2 = obj.UpWind{2}.x * (obj.Mob(:,2) .* rho(:,2));
            dMupyPh2 = obj.UpWind{2}.y * (obj.Mob(:,2) .* rho(:,2));
            dMupzPh2 = obj.UpWind{2}.z * (obj.Mob(:,2) .* rho(:,2));
            vecX1 = min(reshape(obj.U{2}.x(1:Nx,:,:),N,1), 0) .* dMupxPh2; 
            vecX2 = max(reshape(obj.U{2}.x(2:Nx+1,:,:),N,1), 0) .* dMupxPh2; 
            vecY1 = min(reshape(obj.U{2}.y(:,1:Ny,:),N,1), 0) .* dMupyPh2; 
            vecY2 = max(reshape(obj.U{2}.y(:,2:Ny+1,:),N,1), 0) .* dMupyPh2;
            vecZ1 = min(reshape(obj.U{2}.z(:,:,1:Nz),N,1), 0) .* dMupzPh2; 
            vecZ2 = max(reshape(obj.U{2}.z(:,:,2:Nz+1),N,1), 0) .* dMupzPh2;
            acc = pv/dt .* (S(:,2) .* rho(:,2));
            
            DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1,vecZ1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            J1x1ph2 = spdiags(DiagVecs,DiagIndx, N, N);
            
            %% 5. Component 2 x1ph1 block
            J2x1ph1 = -J1x1ph1; 
            
            %% 6. Component 2 x1ph2  block
            J2x1ph2 = -J1x1ph2;
            
            %% 7. Add wells to each block
            [Jp{1}, Jp{2}, JS{1}, JS{2}, ...
                J1x1ph1, J1x1ph2, J2x1ph1, J2x1ph2] = ...
                obj.AddWellsToJacobian(Jp{1}, Jp{2}, JS{1}, JS{2}, J1x1ph1, J1x1ph2, J2x1ph1, J2x1ph2,...
                P, x, rho, ProductionSystem.Wells, ProductionSystem.Reservoir.K);
            
            x1 = x(:,1:2);
            x2 = x(:,3:4);
            
            %% 8. Equilibrium of component 1
            Jeq1p = - spdiags(obj.dKdp(:,1) .* x1(:,2), 0, N, N);
            Jeq1S = speye(N) .* 0;
            Jeq1_x1ph1 = speye(N);
            Jeq1_x1ph2 = - spdiags(obj.K(:,1), 0, N, N);            
                       
            %% 9. Equilibrium of component 2
            Jeq2p = - spdiags(obj.dKdp(:,2) .* x2(:,2), 0, N, N);
            Jeq2S = speye(N) .* 0;
            Jeq2_x1ph1 =  - speye(N);
            Jeq2_x1ph2 = spdiags(obj.K(:,2), 0, N, N);
            
            
            %% 10. Single Phase cells: write dummy equations            
            % Modify transport equation to only solve for dx
            JS{2}(obj.SinglePhase > 0, :) = 0;
            J2x1ph2(obj.SinglePhase == 1,:) = 0;
            J2x1ph1(obj.SinglePhase == 2,:) = 0;
            % Force delta S to be equal to zero
            cells = find(obj.SinglePhase > 0);
            Jeq1S(sub2ind(size(Jeq1S), cells, cells)) = 1;
            Jeq1p (obj.SinglePhase > 0, :) = 0;
            Jeq1_x1ph1 (obj.SinglePhase > 0, :) = 0;
            Jeq1_x1ph2 (obj.SinglePhase > 0, :) = 0;   
            % Force the other dx to be equal to zero
            Jeq2p (obj.SinglePhase > 0, :) = 0;
            % if only vapor force dx1ph2 to be zero
            Jeq2_x1ph1(obj.SinglePhase == 1, :) = 0;
            cells = find(obj.SinglePhase == 1);
            Jeq2_x1ph2(sub2ind(size(Jeq2_x1ph2), cells, cells)) = 1;
            % if only liquid force dx1ph1 to be zero
            Jeq2_x1ph2(obj.SinglePhase == 2, :) = 0;
            cells = find(obj.SinglePhase == 2);
            Jeq2_x1ph1(sub2ind(size(Jeq2_x1ph1), cells, cells)) = 1;
            
            %% Full Jacobian
            Jacobian = [ Jp{1}, JS{1},  J1x1ph1, J1x1ph2;...
                         Jp{2}, JS{2},  J2x1ph1, J2x1ph2;...
                         Jeq1p, Jeq1S, Jeq1_x1ph1, Jeq1_x1ph2;...
                         Jeq2p, Jeq2S, Jeq2_x1ph1, Jeq2_x1ph2];
            
        end
        function [J1p, J2p, J1S, J2S, J1x1ph1, J1x1ph2, J2x1ph1, J2x1ph2] =...
                AddWellsToJacobian(obj, J1p, J2p, J1S, J2S, J1x1ph1, J1x1ph2, J2x1ph1, J2x1ph2, p, x, rho, Wells, K)
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %Injectors
            for i=1:Wells.NofInj
                a = Inj(i).Cells;
                    for j=1:length(a)
                        J1p(a(j),a(j)) = J1p(a(j),a(j)) ...
                            + Inj(i).PI * K(a(j)) * (Inj(i).Mob(:,1) * Inj(i).rho(j,1) * Inj(i).x(:,1) ...
                            + Inj(i).Mob(:,2) *  Inj(i).rho(j,2) * Inj(i).x(:,2));
                        J2p(a(j),a(j)) = J2p(a(j),a(j)) ...
                            + Inj(i).PI * K(a(j)) * (Inj(i).Mob(:,1) * Inj(i).rho(j,1) * Inj(i).x(:,3) ...
                            + Inj(i).Mob(:,2) *  Inj(i).rho(j,2) * Inj(i).x(:,4));
                    end
            end
            
            %Producers
            for i=1:Wells.NofProd
                b = Prod(i).Cells;
                for j=1:length(b)
                    %Pressure blocks
                    J1p(b(j),b(j)) = J1p(b(j),b(j))...
                        + Prod(i).PI * K(b(j)) * (obj.Mob(b(j), 1) * rho(b(j),1) * x(b(j),1) ...
                        + obj.Mob(b(j), 2) * rho(b(j),2) * x(b(j),2))...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * x(b(j),1) * obj.drhodp(b(j),1) * (Prod(i).p(j) - p(b(j))) ...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * x(b(j),2) * obj.drhodp(b(j),2) * (Prod(i).p(j) - p(b(j)));
                    J2p(b(j),b(j)) = J2p(b(j),b(j)) ...
                        + Prod(i).PI * K(b(j)) * (obj.Mob(b(j), 1) * rho(b(j),1) * x(b(j),3) ...
                        + obj.Mob(b(j), 2) * rho(b(j),2) * x(b(j),4))...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * x(b(j),3) * obj.drhodp(b(j),1) * (Prod(i).p(j) - p(b(j))) ...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * x(b(j),4) * obj.drhodp(b(j),2) * (Prod(i).p(j) - p(b(j)));
                    %Saturation blocks
                    J1S(b(j),b(j)) = J1S(b(j),b(j))...
                        - Prod(i).PI * K(b(j)) * (obj.dMob(b(j),1) * rho(b(j),1) * x(b(j),1) ...
                        + obj.dMob(b(j), 2) * rho(b(j), 2) * x(b(j),2)) * (Prod(i).p(j) - p(b(j)));
                    J2S(b(j),b(j)) = J2S(b(j),b(j)) ...
                        - Prod(i).PI * K(b(j)) * (obj.dMob(b(j), 1) * rho(b(j),1) * x(b(j), 3) ...
                        + obj.dMob(b(j),2) * rho(b(j),2) * x(b(j),4)) * (Prod(i).p(j) - p(b(j)));
                    % Mole fractions blocks
                    J1x1ph1(b(j), b(j)) = J1x1ph1(b(j), b(j)) + Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * (Prod(i).p(j) - p(b(j))); 
                    J1x1ph2(b(j), b(j)) = J1x1ph2(b(j), b(j)) + Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * (Prod(i).p(j) - p(b(j)));
                    J2x1ph1(b(j), b(j)) = J2x1ph1(b(j), b(j)) - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * (Prod(i).p(j) - p(b(j)));
                    J2x1ph2(b(j), b(j)) = J2x1ph2(b(j), b(j)) - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * (Prod(i).p(j) - p(b(j)));
                end
                
            end
        end
        function delta = UpdateState(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)   
            %if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                %return
            %else
                Nm =  DiscretizationModel.ReservoirGrid.N;
                %% 1. Update matrix
                % Update Pressure
                Pm = ProductionSystem.Reservoir.State.Properties(['P_', num2str(obj.NofPhases)]);
                Pm.update(delta(1:Nm));
                DeltaLast = zeros(Nm, 1);
                for ph = 1:obj.NofPhases-1
                    Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(ph)]);
                    Sm.update(delta(ph*Nm + 1:(ph+1)*Nm));
                    % Update x unless phase state resulted to be wrong
                    delta1 = delta(obj.NofPhases*Nm + 1:(obj.NofPhases+1)*Nm);
                    delta2 = delta((obj.NofPhases+1)*Nm + 1:(obj.NofPhases+2)*Nm);
                    delta1((Sm.Value > 1)) = 0;
                    delta2((Sm.Value > 1)) = 0;
                    delta1((Sm.Value < 0)) = 0;
                    delta2((Sm.Value < 0)) = 0;
                    xm = ProductionSystem.Reservoir.State.Properties('x_1ph1');
                    xm.update(delta1);
                    xm = ProductionSystem.Reservoir.State.Properties('x_1ph2');
                    xm.update(delta2);
                    xm = ProductionSystem.Reservoir.State.Properties('x_2ph1');
                    xm.update(-delta1);
                    xm = ProductionSystem.Reservoir.State.Properties('x_2ph2');
                    xm.update(-delta2);
            
                    % Single phase from previous solution
                    obj.PreviousSinglePhase = obj.SinglePhase;
                    obj.SinglePhase(Sm.Value > 1) = 1;
                    obj.SinglePhase(Sm.Value < 0) = 2;
                    % Remove values that are not physical
                    Sm.Value = max(Sm.Value, 0);
                    Sm.Value = min(Sm.Value, 1);
                    DeltaLast = DeltaLast + delta(ph*Nm + 1:(ph+1)*Nm);
                end
                Sm = ProductionSystem.Reservoir.State.Properties(['S_', num2str(obj.NofPhases)]);
                Sm.update(-DeltaLast);
                % Remove values that are not physical
                Sm.Value = max(Sm.Value, 0);
                Sm.Value = min(Sm.Value, 1);
                % Update Phase Densities
                FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State, obj.SinglePhase);
                % Update total density
                FluidModel.ComputeTotalDensity(ProductionSystem.Reservoir.State);
                % Compute total mole fractions
                FluidModel.ComputeTotalFractions(ProductionSystem.Reservoir.State, Nm);
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
                        FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(i).State, obj.SinglePhase);
                        % Update total density
                        FluidModel.ComputeTotalDensity(ProductionSystem.FracturesNetwork.Fractures(i).State);
                        % Update Pc
                        FluidModel.ComputePc(ProductionSystem.FracturesNetwork.Fractures(i).State);
                    end
                end
            %end
        end
        function UpdatePandS(obj, delta, Status)
            % Update Solution
            Status.p = Status.p + delta(1:end/4);
            Status.S = Status.S + delta(end/4+1:end/2);
            Status.x1(:,1) = Status.x1(:,1) + delta(end/2 +1:3*end/4);
            Status.x1(:,2) = Status.x1(:,2) + delta(3*end/4 +1:end);
            
            % Single phase from previous solution
            obj.PreviousSinglePhase = obj.SinglePhase;
            obj.SinglePhase(Status.S > 1) = 1;
            obj.SinglePhase(Status.S < 0) = 2;
            
            Status.S = min(Status.S, 1);
            Status.S = max(Status.S, 0);
        end
        function UpdateCompositions(obj, Status, FluidModel, N)
            x(:,1:2) = Status.x1;
            x(:,3:4) = 1 - x(:,1:2);
            
            % THERMODYNAMIC EQUILIBRIUM EQUATIONS
            Req = zeros(N*obj.NofComponents, 1);
            for i=1:obj.NofComponents
                Rcomp = x(:,(i-1)*2 + 1) - obj.K(:,i).*x(:,(i-1)*2 + 2);
                % If Single phase set residual to be zero
                Rcomp(obj.PreviousSinglePhase > 0) = 0;
                Req((i-1)*N+1:i*N) = Rcomp;            
            end
            
            %% 8. Equilibrium of component 1
            Jeq1_x1ph1 = speye(N);
            Jeq1_x1ph2 = - spdiags(obj.K(:,1), 0, N, N);            
                       
            %% 9. Equilibrium of component 2
            Jeq2_x1ph1 =  - speye(N);
            Jeq2_x1ph2 = spdiags(obj.K(:,2), 0, N, N);
            
            %% matrix
            A = [Jeq1_x1ph1, Jeq1_x1ph2;...
                        Jeq2_x1ph1, Jeq2_x1ph2];
            
            delta = A\Req;
            delta1 = delta(1:N);
            %delta1(obj.SinglePhase == 2) = 0;
            delta2 = delta(N+1:end);
            %delta2(obj.SinglePhase == 1) = 0;
            
            Status.x1(:, 1) = Status.x1(:, 1) + delta1;
            Status.x1(:, 2) = Status.x1(:, 2) + delta2;
            
            % Update density
            FluidModel.ComputePhaseDensities(Status);
            
            % Update z
            Status.z = FluidModel.ComputeTotalFractions(Status.S, Status.x1, Status.rho);
            
            % Update Pc
            Status.pc = FluidModel.ComputePc(Status.S);
        end
    end
end