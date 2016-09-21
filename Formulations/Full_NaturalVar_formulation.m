%  Full Natural variable Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 September 2016
%Last modified: 19 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Full_NaturalVar_formulation < fim_formulation
    properties
        NofComponents
        K
        dKdp
        InitialPhaseState
        PreviousSinglePhase
        SinglePhase
    end
    methods
        function obj = Full_NaturalVar_formulation(n_cells, n_components)
            obj.PreviousSinglePhase = zeros(n_cells, 1);
            obj.SinglePhase = zeros(n_cells, 1);
            obj.NofComponents = n_components;
            obj.Tph1 = cell(n_components, 1);
            obj.Tph2 = cell(n_components, 1);
        end
        function SavePhaseState(obj)
            obj.InitialPhaseState = obj.SinglePhase;
        end
        function Reset(obj)
            obj.SinglePhase = obj.InitialPhaseState;
            obj.PreviousSinglePhase = obj.InitialPhaseState;
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.S);
            obj.dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State.S);
            obj.drho = FluidModel.DrhoDp(ProductionSystem.Reservoir.State.p);
            obj.dPc = FluidModel.DPcDS(ProductionSystem.Reservoir.State.S);
            obj.K = FluidModel.KvaluesCalculator.Compute(ProductionSystem.Reservoir.State.p, ProductionSystem.Reservoir.State.T, FluidModel.Components);
            obj.dKdp = FluidModel.KvaluesCalculator.DKvalDp(ProductionSystem.Reservoir.State.p);
            obj.SinglePhase = FluidModel.CheckNumberOfPhases(obj.SinglePhase, obj.PreviousSinglePhase, ProductionSystem.Reservoir.State.z, obj.K);
         end
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)                   
            %Create local variables
            N = DiscretizationModel.ReservoirGrid.N;
            s_old = State0.S;
            x_old(:,1:2) = State0.x1;
            x_old(:,3:4) = 1 - x_old(:,1:2);
            s = ProductionSystem.Reservoir.State.S;
            x(:,1:2) = ProductionSystem.Reservoir.State.x1;
            x(:,3:4) = 1 - x(:,1:2);
            s2 = 1 - s;
            s2_old = 1 - s_old;
            % Phase potentials
            Pot_ph1 = ProductionSystem.Reservoir.State.Pot(:,1);
            Pot_ph2 = ProductionSystem.Reservoir.State.Pot(:,2);
            % Pore volume
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            %Density
            rho = ProductionSystem.Reservoir.State.rho;
            rho_old = State0.rho;
            
            %Accumulation term
            A = speye(N)*pv/dt;
                     
            %Source terms
            q = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            
            %% RESIDUAL
            Rbalance = zeros(N*obj.NofComponents, 1);
            Req = zeros(N*obj.NofComponents, 1);
            for i=1:obj.NofComponents
                % 1. MASS CONSERVATION EQUATIONS                 
                m = x(:,(i-1)*2+1) .* rho(:,1) .* s + x(:,(i-1)*2+2) .* rho(:,2) .* s2;
                m_old = x_old(:,(i-1)*2+1) .* rho_old(:,1) .* s_old + x_old(:,(i-1)*2+2) .* rho_old(:,2) .* s2_old;
                % Phase Transmissibilities
                obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, rho, x(:,(i-1)*2+1:(i-1)*2+2), i);          
                % Residual
                Rbalance((i-1)*N+1:i*N) = ...
                    A * m - A * m_old...          %Accumulation term
                    + obj.Tph1{i} *  Pot_ph1 ...     %Convective term
                    + obj.Tph2{i} *  Pot_ph2...
                    - q(:,i);                           %Source terms
                
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
            N = DiscretizationModel.ReservoirGrid.N;
            pv = DiscretizationModel.ReservoirGrid.Volume*ProductionSystem.Reservoir.Por;
            x(:,1:2) = ProductionSystem.Reservoir.State.x1;
            x(:,3:4) = 1 - x(:,1:2);
            x1 = ProductionSystem.Reservoir.State.x1;
            x2 = 1 - x1;
            s = ProductionSystem.Reservoir.State.S;
            rho = ProductionSystem.Reservoir.State.rho;
            
            % Fill in block by block
            Jp = cell(obj.NofComponents, 1);
            JS = cell(obj.NofComponents, 1);
            for i=1:obj.NofComponents
                %% 1. Component i pressure block
                % 1.a: divergence
                Jp{i} = obj.Tph1{i}  + obj.Tph2{i};
                % 1.b: compressibility part
                dMupxPh1 = obj.UpWindPh1.x * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drho(:,1));
                dMupyPh1 = obj.UpWindPh1.y * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drho(:,1));
                dMupxPh2 = obj.UpWindPh2.x * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drho(:,2));
                dMupyPh2 = obj.UpWindPh2.y * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drho(:,2));
                
                vecX1 = min(reshape(obj.Uph1.x(1:Nx,:),N,1), 0).*dMupxPh1 + min(reshape(obj.Uph2.x(1:Nx,:),N,1), 0).*dMupxPh2;
                vecX2 = max(reshape(obj.Uph1.x(2:Nx+1,:),N,1), 0).*dMupxPh1 + max(reshape(obj.Uph2.x(2:Nx+1,:),N,1), 0).*dMupxPh2;
                vecY1 = min(reshape(obj.Uph1.y(:,1:Ny),N,1), 0).*dMupyPh1 + min(reshape(obj.Uph2.y(:,1:Ny),N,1), 0).*dMupyPh2;
                vecY2 = max(reshape(obj.Uph1.y(:,2:Ny+1),N,1), 0).*dMupyPh1 + max(reshape(obj.Uph2.y(:,2:Ny+1),N,1), 0).*dMupyPh2;
                acc = pv/dt .* ( x(:,(i-1)*2+1) .* obj.drho(:,1) .* s + x(:,(i-1)*2+2) .* obj.drho(:,2) .* (1-s));
                DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
                DiagIndx = [-Nx, -1, 0, 1, Nx];
                Jp{i} = Jp{i} + spdiags(DiagVecs, DiagIndx, N, N);
                
                %% 2. Component i saturation block
                dMupxPh1 = obj.UpWindPh1.x*(obj.dMob(:,1) .* x(:,(i-1)*2+1) .* rho(:,1));
                dMupyPh1 = obj.UpWindPh1.y*(obj.dMob(:,1) .* x(:,(i-1)*2+1) .* rho(:,1));
                dMupxPh2 = obj.UpWindPh2.x*(obj.dMob(:,2) .* x(:,(i-1)*2+2) .* rho(:,2));
                dMupyPh2 = obj.UpWindPh2.y*(obj.dMob(:,2) .* x(:,(i-1)*2+2) .* rho(:,2));
                vecX1 = min(reshape(obj.Uph1.x(1:Nx,:),N,1), 0).*dMupxPh1 + min(reshape(obj.Uph2.x(1:Nx,:),N,1), 0).*dMupxPh2;
                vecX2 = max(reshape(obj.Uph1.x(2:Nx+1,:),N,1), 0).*dMupxPh1 + max(reshape(obj.Uph2.x(2:Nx+1,:),N,1), 0).*dMupxPh2;
                vecY1 = min(reshape(obj.Uph1.y(:,1:Ny),N,1), 0).*dMupyPh1 + min(reshape(obj.Uph2.y(:,1:Ny),N,1), 0).*dMupyPh2;
                vecY2 = max(reshape(obj.Uph1.y(:,2:Ny+1),N,1), 0).*dMupyPh1 + max(reshape(obj.Uph2.y(:,2:Ny+1),N,1), 0).*dMupyPh2;
                acc = pv/dt .* (x(:,(i-1)*2+1) .* rho(:,1) - x(:,(i-1)*2+2) .* rho(:,2));
                DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
                DiagIndx = [-Nx, -1, 0, 1, Nx];
                JS{i} = spdiags(DiagVecs,DiagIndx, N, N);
                % Capillarity
                JS{i} = JS{i} - obj.Tph1{i} * spdiags(obj.dPc, 0, N, N);               
            end
            
            %% 3. Component 1 x1ph1 block
            dMupxPh1 = obj.UpWindPh1.x * (obj.Mob(:,1) .* rho(:,1));
            dMupyPh1 = obj.UpWindPh1.y * (obj.Mob(:,1) .* rho(:,1));
            vecX1 = min(reshape(obj.Uph1.x(1:Nx,:),N,1), 0).*dMupxPh1; 
            vecX2 = max(reshape(obj.Uph1.x(2:Nx+1,:),N,1), 0).*dMupxPh1; 
            vecY1 = min(reshape(obj.Uph1.y(:,1:Ny),N,1), 0).*dMupyPh1; 
            vecY2 = max(reshape(obj.Uph1.y(:,2:Ny+1),N,1), 0).*dMupyPh1; 
            acc = pv/dt .* (s .* rho(:,1));
            DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
            DiagIndx = [-Nx, -1, 0, 1, Nx];
            J1x1ph1 = spdiags(DiagVecs,DiagIndx, N, N);
            
            %% 4. Component 1 x1ph2  block
            dMupxPh2 = obj.UpWindPh1.x * (obj.Mob(:,2) .* rho(:,2));
            dMupyPh2 = obj.UpWindPh1.y * (obj.Mob(:,2) .* rho(:,2));
            vecX1 = min(reshape(obj.Uph2.x(1:Nx,:),N,1), 0) .* dMupxPh2; 
            vecX2 = max(reshape(obj.Uph2.x(2:Nx+1,:),N,1), 0) .* dMupxPh2; 
            vecY1 = min(reshape(obj.Uph2.y(:,1:Ny),N,1), 0) .* dMupyPh2; 
            vecY2 = max(reshape(obj.Uph2.y(:,2:Ny+1),N,1), 0) .* dMupyPh2; 
            acc = pv/dt .* ((1 - s) .* rho(:,2));
            DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
            DiagIndx = [-Nx, -1, 0, 1, Nx];
            J1x1ph2 = spdiags(DiagVecs,DiagIndx, N, N);
            
            %% 5. Component 2 x1ph1 block
            J2x1ph1 = -J1x1ph1; 
            
            %% 6. Component 2 x1ph2  block
            J2x1ph2 = -J1x1ph2;
            
            %% 7. Add wells to each block
            [Jp{1}, Jp{2}, JS{1}, JS{2}, ...
                J1x1ph1, J1x1ph2, J2x1ph1, J2x1ph2] = ...
                obj.AddWellsToJacobian(Jp{1}, Jp{2}, JS{1}, JS{2}, J1x1ph1, J1x1ph2, J2x1ph1, J2x1ph2,...
                ProductionSystem.Reservoir.State, ProductionSystem.Wells, ProductionSystem.Reservoir.K);
            
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
        function delta = UpdateState(obj, delta, Status, FluidModel)
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
            
            % Update density
            for i=1:FluidModel.NofPhases
                Status.rho(:, i) = FluidModel.Phases(i).ComputeDensity(Status.p);
            end
            
            % Update z
            Status.z = FluidModel.ComputeTotalFractions(Status.S, Status.x1, Status.rho);
            
            % Update Pc
            Status.pc = FluidModel.ComputePc(Status.S);
            
            % Compute phase potential
            h = FluidModel.GravityModel.ComputeBuoyancyTerm(Status.rho);
            Status.ComputePhasePotential(h);
        end
        function  TransmissibilityMatrix(obj, Grid, Rho, x, i)
            %%%Transmissibility matrix construction
            Tx = zeros(Grid.Nx+1, Grid.Ny);
            Ty = zeros(Grid.Nx, Grid.Ny+1);
            
            %% PHASE 1 
            %Apply upwind operator
            Mupx = obj.UpWindPh1.x*(obj.Mob(:,1) .* Rho(:,1) .* x(:,1)); 
            Mupy = obj.UpWindPh1.y*(obj.Mob(:,1) .* Rho(:,1) .* x(:,1));
            Mupx = reshape(Mupx, Grid.Nx, Grid.Ny);
            Mupy = reshape(Mupy, Grid.Nx, Grid.Ny);
            Tx(2:Grid.Nx,:)= Grid.Tx(2:Grid.Nx,:).*Mupx(1:Grid.Nx-1,:);
            Ty(:,2:Grid.Ny)= Grid.Ty(:,2:Grid.Ny).*Mupy(:,1:Grid.Ny-1);
            
            %Construct matrix
            x1 = reshape(Tx(1:Grid.Nx,:), Grid.N, 1);
            x2 = reshape(Tx(2:Grid.Nx+1,:), Grid.N, 1);
            y1 = reshape(Ty(:,1:Grid.Ny), Grid.N, 1);
            y2 = reshape(Ty(:,2:Grid.Ny+1), Grid.N, 1);
            DiagVecs = [-y2,-x2,y2+x2+y1+x1,-x1,-y1];
            DiagIndx = [-Grid.Nx, -1, 0, 1, Grid.Nx];
            obj.Tph1{i} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
            
            %% PHASE 2 
            %Apply upwind operator
            Mupx = obj.UpWindPh2.x*(obj.Mob(:,2) .* Rho(:,2) .* x(:,2));
            Mupy = obj.UpWindPh2.y*(obj.Mob(:,2) .* Rho(:,2) .* x(:,2));
            Mupx = reshape(Mupx, Grid.Nx, Grid.Ny);
            Mupy = reshape(Mupy, Grid.Nx, Grid.Ny);
            Tx(2:Grid.Nx,:)= Grid.Tx(2:Grid.Nx,:).*Mupx(1:Grid.Nx-1,:);
            Ty(:,2:Grid.Ny)= Grid.Ty(:,2:Grid.Ny).*Mupy(:,1:Grid.Ny-1);
            
            %Construct matrix
            x1 = reshape(Tx(1:Grid.Nx,:), Grid.N, 1);
            x2 = reshape(Tx(2:Grid.Nx+1,:), Grid.N, 1);
            y1 = reshape(Ty(:,1:Grid.Ny), Grid.N, 1);
            y2 = reshape(Ty(:,2:Grid.Ny+1), Grid.N, 1);
            DiagVecs = [-y2,-x2,y2+x2+y1+x1,-x1,-y1];
            DiagIndx = [-Grid.Nx, -1, 0, 1, Grid.Nx];
            obj.Tph2{i} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
        end
        function q = ComputeSourceTerms(obj, N, Wells)
            q = zeros(N, 2);
            %Injectors
            for i=1:Wells.NofInj
                c = Wells.Inj(i).Cells;
                q(c, :) = Wells.Inj(i).QComponents(:,:);
            end
            %Producers
            for i=1:Wells.NofProd
                c = Wells.Prod(i).Cells;
                q(c, :) = Wells.Prod(i).QComponents(:,:);
            end
        end
        function [J1p, J2p, J1S, J2S, J1x1ph1, J1x1ph2, J2x1ph1, J2x1ph2] =...
                AddWellsToJacobian(obj, J1p, J2p, J1S, J2S, J1x1ph1, J1x1ph2, J2x1ph1, J2x1ph2, State, Wells, K)
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %Injectors
            for i=1:Wells.NofInj
                a = Inj(i).Cells;
                    for j=1:length(a)
                        J1p(a(j),a(j)) = J1p(a(j),a(j)) ...
                            + Inj(i).PI * K(a(j)) * (Inj(i).Mob(:,1) * Inj(i).rho(:,1) * Inj(i).x1(:,1) ...
                            + Inj(i).Mob(:,2) *  Inj(i).rho(:,2) * Inj(i).x1(:,2));
                        J2p(a(j),a(j)) = J2p(a(j),a(j)) ...
                            + Inj(i).PI * K(a(j)) * (Inj(i).Mob(:,1) * Inj(i).rho(:,1) * Inj(i).x2(:,1) ...
                            + Inj(i).Mob(:,2) *  Inj(i).rho(:,2) * Inj(i).x2(:,2));
                    end
            end
            
            %Producers
            for i=1:Wells.NofProd
                b = Prod(i).Cells;
                for j=1:length(b)
                    %Pressure blocks
                    J1p(b(j),b(j)) = J1p(b(j),b(j))...
                        + Prod(i).PI * K(b(j)) * (obj.Mob(b(j), 1) * State.rho(b(j),1) * State.x1(b(j),1) ...
                        + obj.Mob(b(j), 2) * State.rho(b(j),2) * State.x1(b(j),2))...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * State.x1(b(j),1) * obj.drho(b(j),1) * (Prod(i).p - State.p(b(j))) ...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * State.x1(b(j),2) * obj.drho(b(j),2) * (Prod(i).p - State.p(b(j)));
                    J2p(b(j),b(j)) = J2p(b(j),b(j)) ...
                        + Prod(i).PI * K(b(j)) * (obj.Mob(b(j), 1) * State.rho(b(j),1) * (1 - State.x1(b(j),1)) ...
                        + obj.Mob(b(j), 2) * State.rho(b(j),2) * (1 - State.x1(b(j),2)))...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * (1 - State.x1(b(j),1)) * obj.drho(b(j),1) * (Prod(i).p - State.p(b(j))) ...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * (1 - State.x1(b(j),2)) * obj.drho(b(j),2) * (Prod(i).p - State.p(b(j)));
                    %Saturation blocks
                    J1S(b(j),b(j)) = J1S(b(j),b(j))...
                        - Prod(i).PI * K(b(j)) * (obj.dMob(b(j),1) * State.rho(b(j),1) * State.x1(b(j),1) ...
                        + obj.dMob(b(j), 2) * State.rho(b(j), 2) * State.x1(b(j),2)) * (Prod(i).p - State.p(b(j)));
                    J2S(b(j),b(j)) = J2S(b(j),b(j)) ...
                        - Prod(i).PI * K(b(j)) * (obj.dMob(b(j), 1) * State.rho(b(j),1) * (1- State.x1(b(j), 1)) ...
                        + obj.dMob(b(j),2) * State.rho(b(j),2) * (1 - State.x1(b(j),2))) * (Prod(i).p - State.p(b(j)));
                    % Mole fractions blocks
                    J1x1ph1(b(j), b(j)) = J1x1ph1(b(j), b(j)) + Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * (Prod(i).p - State.p(b(j))); 
                    J1x1ph2(b(j), b(j)) = J1x1ph2(b(j), b(j)) + Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * (Prod(i).p - State.p(b(j)));
                    J2x1ph1(b(j), b(j)) = J2x1ph1(b(j), b(j)) - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * (Prod(i).p - State.p(b(j)));
                    J2x1ph2(b(j), b(j)) = J2x1ph2(b(j), b(j)) - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * (Prod(i).p - State.p(b(j)));
                    
                end
                
            end
        end
    end
end