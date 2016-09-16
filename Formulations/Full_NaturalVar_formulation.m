%  Natural variable Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Full_NaturalVar_formulation < fim_formulation
    properties
        Tc1
        Tc2
        Tph1
        Tph2
        Mob
        dMob
        dPc
        drho
        K
        dKdp
        SinglePhase = zeros(5, 1);
    end
    methods
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.S);
            obj.dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State.S);
            obj.drho = FluidModel.DrhoDp(ProductionSystem.Reservoir.State.p);
            obj.dPc = FluidModel.DPcDS(ProductionSystem.Reservoir.State.S);
            obj.K = FluidModel.KvaluesCalculator.Compute(ProductionSystem.Reservoir.State.p, ProductionSystem.Reservoir.State.T, FluidModel.Components);
            obj.dKdp = FluidModel.KvaluesCalculator.DKvalDp(ProductionSystem.Reservoir.State.p);
            
         end
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)                   
            %Create local variables
            N = DiscretizationModel.ReservoirGrid.N;
            s_old = State0.S;
            x1_old = State0.x1;
            x2_old = 1 - x1_old;
            p = ProductionSystem.Reservoir.State.p;
            s = ProductionSystem.Reservoir.State.S;
            x1 = ProductionSystem.Reservoir.State.x1;
            x2 = 1 - x1;
            s2 = 1 - s;
            s2_old = 1 - s_old;
            Pc = ProductionSystem.Reservoir.State.pc;
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            %Density
            rho = ProductionSystem.Reservoir.State.rho;
            rho_old = State0.rho;
            
            % 1. MASS CONSERVATION EQUATIONS
            
            %Accumulation term
            A = speye(N)*pv/dt;
            %Component 1
            m1 = x1(:,1) .* rho(:,1) .* s + x1(:,2) .* rho(:,2) .* s2;
            m1_old = x1_old(:,1) .* rho_old(:,1) .* s_old + x1_old(:,2) .* rho_old(:,2) .* s2_old;
            %Component 2
            m2 = x2(:,1) .* rho(:,1) .* s + x2(:,2) .* rho(:,2) .* s2;
            m2_old = x2_old(:,1) .* rho_old(:,1).* s_old + x2_old(:,2) .* rho_old(:,2) .* s2_old;
            
            %Convective term
            obj.Tc1 = obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, rho, x1);
            obj.Tc2 = obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, rho, x2);
            
            %Capillarity
            obj.Tph1 = obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, rho, [ones(N,1), zeros(N,1)]);
            
            %Gravity
            G = obj.ComputeGravityTerm(N);
            
            %Source terms
            [q1, q2] = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            
            %% RESIDUAL
            %Component 1
            R1 = A * m1 - A * m1_old...           %Accumulation term
                + obj.Tc1 * p...                       %Convective term
                - x1(:,1) .* (obj.Tph1*Pc)...           %Capillarity
                + G*p...                          %Gravity
                - q1;                             %Source terms
            %Component 2
            R2 = A * m2 - A * m2_old...
                + obj.Tc2 * p...                       %Convective term
                - x2(:,1) .* (obj.Tph1*Pc)...           %Capillarity
                + G*p...                          %Gravity
                - q2;                             %Source terms
            
            
           % 2. THERMODYNAMIC EQUILIBRIUM EQUATIONS
           Req1 = x1(:,1) - obj.K(:,1).*x1(:,2);
           Req2 = x2(:,1) - obj.K(:,2).*x2(:,2);
           
           % If Single phase set residual to be zero
           R2(obj.SinglePhase > 0) = 0;
           Req1(obj.SinglePhase > 0) = 0;
           Req2(obj.SinglePhase > 0) = 0;
           
            %Stick them together
            Residual = [R1; R2; Req1; Req2];
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            %Build FIM Jacobian
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            N = DiscretizationModel.ReservoirGrid.N;
            pv = DiscretizationModel.ReservoirGrid.Volume*ProductionSystem.Reservoir.Por;
            x1 = ProductionSystem.Reservoir.State.x1;
            x2 = 1 - x1;
            s = ProductionSystem.Reservoir.State.S;
            rho = ProductionSystem.Reservoir.State.rho;
            
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            
            %% 1. Component 1 pressure block            
            % 1.a: divergence
            J1p = obj.Tc1; 
            % 1.b: compressibility part
            dMupxPh1 = obj.UpWindPh1.x * (obj.Mob(:, 1) .* x1(:,1) .* obj.drho(:,1));
            dMupyPh1 = obj.UpWindPh1.y * (obj.Mob(:, 1) .* x1(:,1) .* obj.drho(:,1));
            dMupxPh2 = obj.UpWindPh2.x * (obj.Mob(:, 2) .* x1(:,2) .* obj.drho(:,2));
            dMupyPh2 = obj.UpWindPh2.y * (obj.Mob(:, 2) .* x1(:,2) .* obj.drho(:,2));
            
            vecX1 = min(reshape(obj.Uph1.x(1:Nx,:),N,1), 0).*dMupxPh1 + min(reshape(obj.Uph2.x(1:Nx,:),N,1), 0).*dMupxPh2;
            vecX2 = max(reshape(obj.Uph1.x(2:Nx+1,:),N,1), 0).*dMupxPh1 + max(reshape(obj.Uph2.x(2:Nx+1,:),N,1), 0).*dMupxPh2;
            vecY1 = min(reshape(obj.Uph1.y(:,1:Ny),N,1), 0).*dMupyPh1 + min(reshape(obj.Uph2.y(:,1:Ny),N,1), 0).*dMupyPh2;
            vecY2 = max(reshape(obj.Uph1.y(:,2:Ny+1),N,1), 0).*dMupyPh1 + max(reshape(obj.Uph2.y(:,2:Ny+1),N,1), 0).*dMupyPh2;
            acc = pv/dt .* ( x1(:,1) .* obj.drho(:,1) .* s + x1(:,2) .* obj.drho(:,2) .* (1-s));
            DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
            DiagIndx = [-Nx, -1, 0, 1, Nx];
            J1p = J1p + spdiags(DiagVecs, DiagIndx, N, N);
            
            %% 2. Component 2 pressure block          
            % 2.a divergence
            J2p = obj.Tc2;
            % 2.b: compressibility part
            dMupxPh1 = obj.UpWindPh1.x*(obj.Mob(:, 1) .* x2(:,1) .* obj.drho(:,1));
            dMupyPh1 = obj.UpWindPh1.y*(obj.Mob(:, 1) .* x2(:,1) .* obj.drho(:,1));
            dMupxPh2 = obj.UpWindPh2.x*(obj.Mob(:, 2) .* x2(:,2) .* obj.drho(:,2));
            dMupyPh2 = obj.UpWindPh2.y*(obj.Mob(:, 2) .* x2(:,2) .* obj.drho(:,2));
            
            vecX1 = min(reshape(obj.Uph1.x(1:Nx,:),N,1), 0).*dMupxPh1 + min(reshape(obj.Uph2.x(1:Nx,:),N,1),0).*dMupxPh2;
            vecX2 = max(reshape(obj.Uph1.x(2:Nx+1,:),N,1), 0).*dMupxPh1 + max(reshape(obj.Uph2.x(2:Nx+1,:),N,1),0).*dMupxPh2;
            vecY1 = min(reshape(obj.Uph1.y(:,1:Ny),N,1), 0).*dMupyPh1 + min(reshape(obj.Uph2.y(:,1:Ny),N,1),0).*dMupyPh2;
            vecY2 = max(reshape(obj.Uph1.y(:,2:Ny+1),N,1), 0).*dMupyPh1 + max(reshape(obj.Uph2.y(:,2:Ny+1),N,1),0).*dMupyPh2;
            acc = pv/dt .* (x2(:,1) .* obj.drho(:,1) .* s + x2(:,2) .* obj.drho(:,2) .* (1-s));
            DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
            DiagIndx = [-Nx, -1, 0, 1, Nx];
            J2p = J2p + spdiags(DiagVecs, DiagIndx, N, N);
            
            %% 3. Component 1 saturation block
            dMupxPh1 = obj.UpWindPh1.x*(obj.dMob(:,1) .* x1(:,1) .* rho(:,1));
            dMupyPh1 = obj.UpWindPh1.y*(obj.dMob(:,1) .* x1(:,1) .* rho(:,1));
            dMupxPh2 = obj.UpWindPh2.x*(obj.dMob(:,2) .* x1(:,2) .* rho(:,2));
            dMupyPh2 = obj.UpWindPh2.y*(obj.dMob(:,2) .* x1(:,2) .* rho(:,2));
            vecX1 = min(reshape(obj.Uph1.x(1:Nx,:),N,1), 0).*dMupxPh1 + min(reshape(obj.Uph2.x(1:Nx,:),N,1), 0).*dMupxPh2;
            vecX2 = max(reshape(obj.Uph1.x(2:Nx+1,:),N,1), 0).*dMupxPh1 + max(reshape(obj.Uph2.x(2:Nx+1,:),N,1), 0).*dMupxPh2;
            vecY1 = min(reshape(obj.Uph1.y(:,1:Ny),N,1), 0).*dMupyPh1 + min(reshape(obj.Uph2.y(:,1:Ny),N,1), 0).*dMupyPh2;
            vecY2 = max(reshape(obj.Uph1.y(:,2:Ny+1),N,1), 0).*dMupyPh1 + max(reshape(obj.Uph2.y(:,2:Ny+1),N,1), 0).*dMupyPh2;
            acc = pv/dt .* (x1(:,1) .* rho(:,1) - x1(:,2) .* rho(:,2));
            DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
            DiagIndx = [-Nx, -1, 0, 1, Nx];
            J1S = spdiags(DiagVecs,DiagIndx, N, N);
            % Capillarity
            J1S = J1S - obj.Tph1 * spdiags(x1(:,1).* obj.dPc, 0, N, N);
            
            %% 4. Component 2 saturation block
            dMupxPh1 = obj.UpWindPh1.x*(obj.dMob(:, 1) .* x2(:,1) .* rho(:,1));
            dMupyPh1 = obj.UpWindPh1.y*(obj.dMob(:, 1) .* x2(:,1) .* rho(:,1));
            dMupxPh2 = obj.UpWindPh2.x*(obj.dMob(:, 2) .* x2(:,2) .* rho(:,2));
            dMupyPh2 = obj.UpWindPh2.y*(obj.dMob(:, 2) .* x2(:,2) .* rho(:,2));
            vecX1 = min(reshape(obj.Uph1.x(1:Nx,:),N,1), 0).*dMupxPh1 + min(reshape(obj.Uph2.x(1:Nx,:),N,1), 0).*dMupxPh2;
            vecX2 = max(reshape(obj.Uph1.x(2:Nx+1,:),N,1), 0).*dMupxPh1 + max(reshape(obj.Uph2.x(2:Nx+1,:),N,1), 0).*dMupxPh2;
            vecY1 = min(reshape(obj.Uph1.y(:,1:Ny),N,1), 0).*dMupyPh1 + min(reshape(obj.Uph2.y(:,1:Ny),N,1), 0).*dMupyPh2;
            vecY2 = max(reshape(obj.Uph1.y(:,2:Ny+1),N,1), 0).*dMupyPh1 + max(reshape(obj.Uph2.y(:,2:Ny+1),N,1), 0).*dMupyPh2;
            acc = pv/dt .* (x2(:,1) .* rho(:,1) - x2(:,2) .* rho(:,2));
            DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
            DiagIndx = [-Nx, -1, 0, 1, Nx];
            J2S = spdiags(DiagVecs, DiagIndx, N, N);
            % Capillarity
            J2S = J2S - obj.Tph1 * spdiags(x2(:,1).* obj.dPc, 0, N, N);
            
            
            %% 5. Component 1 x1ph1 block
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
            
            %% 6. Component 1 x1ph2  block
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
            
            %% 7. Component 2 x1ph1 block
            J2x1ph1 = -J1x1ph1; 
            
            %% 8. Component 2 x1ph2  block
            J2x1ph2 = -J1x1ph2;
            
            %% Add wells to each block
            [J1p, J2p, J1S, J2S, ...
                J1x1ph1, J1x1ph2, J2x1ph1, J2x1ph2] = ...
                obj.AddWellsToJacobian(J1p, J2p, J1S, J2S, J1x1ph1, J1x1ph2, J2x1ph1, J2x1ph2,...
                ProductionSystem.Reservoir.State, ProductionSystem.Wells, ProductionSystem.Reservoir.K);
            
            %% 9. Equilibrium of component 1
            Jeq1p = - spdiags(obj.dKdp(:,1) .* x1(:,2), 0, N, N);
            Jeq1S = speye(N) .* 0;
            Jeq1_x1ph1 = speye(N);
            Jeq1_x1ph2 = - spdiags(obj.K(:,1), 0, N, N);            
                       
            %% 10. Equilibrium of component 2
            Jeq2p = - spdiags(obj.dKdp(:,2) .* x2(:,2), 0, N, N);
            Jeq2S = speye(N) .* 0;
            Jeq2_x1ph1 =  - speye(N);
            Jeq2_x1ph2 = spdiags(obj.K(:,2), 0, N, N);
            
            
            %% Single Phase cells: write dummy equations
            
            % Force dS to be equal to 0
            J2p(obj.SinglePhase > 0, :) = 0;
            J2S(obj.SinglePhase > 0, :) = 0;
            J2x1ph1(obj.SinglePhase > 0, :) = 0;
            J2x1ph2(obj.SinglePhase > 0, :) = 0;
            cells = find(obj.SinglePhase > 0);
            J2S(sub2ind(size(J2S), cells, cells)) = 1;

            % Force dxs to be equal to zero for the face that s not there
            Jeq1p(obj.SinglePhase == 2,:) = 0;
            Jeq1_x1ph2(obj.SinglePhase == 2,:) = 0;
            Jeq2p(obj.SinglePhase == 1,:) = 0;
            Jeq2_x1ph1(obj.SinglePhase == 1,:) = 0;
            cells = find(obj.SinglePhase == 1);
            Jeq2_x1ph2(sub2ind(size(Jeq2_x1ph2), cells, cells)) = 1;
            
            %% Full Jacobian
            Jacobian = [ J1p, J1S,  J1x1ph1, J1x1ph2;...
                         J2p, J2S,  J2x1ph1, J2x1ph2;...
                         Jeq1p, Jeq1S, Jeq1_x1ph1, Jeq1_x1ph2;...
                         Jeq2p, Jeq2S, Jeq2_x1ph1, Jeq2_x1ph2];
            
        end
        function delta = UpdateState(obj, delta, Status, FluidModel)
            % Update Solution
            Status.p = Status.p + delta(1:end/4);
            Status.S = Status.S + delta(end/4+1:end/2);
            Status.x1(:,1) = Status.x1(:,1) + delta(end/2 +1:3*end/4);
            Status.x1(:,2) = Status.x1(:,2) + delta(3*end/4 +1:end);
            
            obj.SinglePhase(Status.S > 1) = 1;
            obj.SinglePhase(Status.S < 0) = 2;
            Status.x1(Status.S > 1, 2) = 1;
            Status.x1(Status.S < 0, 1) = 1;
            
            Status.S = min(Status.S, 1);
            Status.S = max(Status.S, 0);
            
            % If Single phase
            Status.S(obj.SinglePhase == 1) = 1;
            Status.S(obj.SinglePhase == 2) = 0;
            
            % Update density
            for i=1:FluidModel.NofPhases
                Status.rho(:, i) = FluidModel.Phases(i).ComputeDensity(Status.p);
            end
            
            % Update z
            Status.z = FluidModel.ComputeMassFractions(Status.S, Status.x1, Status.rho);
        end
        function  T = TransmissibilityMatrix(obj, Grid, Rho, x)
            %%%Transmissibility matrix construction
            Tx = zeros(Grid.Nx+1, Grid.Ny);
            Ty = zeros(Grid.Nx, Grid.Ny+1);
            
            %Apply upwind operator
            Mupx = obj.UpWindPh1.x*(obj.Mob(:,1) .* Rho(:,1) .* x(:,1)) + obj.UpWindPh2.x*(obj.Mob(:,2) .* Rho(:,2) .* x(:,2));
            Mupy = obj.UpWindPh1.y*(obj.Mob(:,1) .* Rho(:,1) .* x(:,1)) + obj.UpWindPh2.y*(obj.Mob(:,2) .* Rho(:,2) .* x(:,2));
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
            T = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
        end
        function [q1, q2] = ComputeSourceTerms(obj, N, Wells)
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
            q1 = q(:,1);
            q2 = q(:,2);
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
        function G = ComputeGravityTerm(obj, N)
               G = zeros(N);
        end
    end
end