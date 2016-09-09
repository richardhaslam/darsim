%  Natural variable Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef NaturalVar_formulation < fim_formulation
    properties
        T1
        T2
        Tw
        Tnw
        Mob
        dMob
        dPc
        drho
    end
    methods
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.S);
            obj.dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State.S);
            obj.drho = FluidModel.DrhoDp(ProductionSystem.Reservoir.State.p);
            %obj.dPc = FluidModel.CapillaryModel.Derivative(ProductionSystem.Reservoir.State.S);
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
            
            %Accumulation term
            A = speye(N)*pv/dt;
            %Component 1
            m1 = x1(:,1) .* rho(:,1) .* s + x1(:,2) .* rho(:,2) .* s2;
            m1_old = x1_old(:,1) .* rho_old(:,1) .* s_old + x1_old(:,2) .* rho_old(:,2) .* s2_old;
            %Component 2
            m2 = x2(:,1) .* rho(:,1) .* s + x2(:,2) .* rho(:,2) .* s2;
            m2_old = x2_old(:,1) .* rho_old(:,1).* s_old + x2_old(:,2) .* rho_old(:,2) .* s2_old;
            
            %Convective term
            obj.T1 = obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, rho, x1);
            obj.T2 = obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, rho, x2);
            
            %Capillarity
            obj.Tw = obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, rho, [ones(N,1), zeros(N,1)]);
            
            %Gravity
            G = obj.ComputeGravityTerm(N);
            
            %Source terms
            [q1, q2] = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            
            %% RESIDUAL
            %Component 1
            R1 = A * m1 - A * m1_old...           %Accumulation term
                + obj.T1 * p...                       %Convective term
                - x1(:,1) .* (obj.Tw*Pc)...           %Capillarity
                + G*p...                          %Gravity
                - q1;                             %Source terms
            %Component 2
            R2 = A * m2 - A * m2_old...
                + obj.T2 * p...                       %Convective term
                - x2(:,1) .* (obj.Tw*Pc)...           %Capillarity
                + G*p...                          %Gravity
                - q2;                             %Source terms
            
            %Stick them together
            Residual = [R1; R2];
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
            
            %% 1 Component 1 pressure block
            
            % 1.a: divergence
            J1p = obj.T1;
            
            % 1.b: compressibility part
            dMupxW = obj.UpWindW.x * (obj.Mob(:, 1) .* x1(:,1) .* obj.drho(:,1));
            dMupyW = obj.UpWindW.y * (obj.Mob(:, 1) .* x1(:,1) .* obj.drho(:,1));
            dMupxNw = obj.UpWindNw.x * (obj.Mob(:, 2) .* x1(:,2) .* obj.drho(:,2));
            dMupyNw = obj.UpWindNw.y * (obj.Mob(:, 2) .* x1(:,2) .* obj.drho(:,2));
            
            vecX1 = min(reshape(obj.Uw.x(1:Nx,:),N,1), 0).*dMupxW + min(reshape(obj.Unw.x(1:Nx,:),N,1), 0).*dMupxNw;
            vecX2 = max(reshape(obj.Uw.x(2:Nx+1,:),N,1), 0).*dMupxW + max(reshape(obj.Unw.x(2:Nx+1,:),N,1), 0).*dMupxNw;
            vecY1 = min(reshape(obj.Uw.y(:,1:Ny),N,1), 0).*dMupyW + min(reshape(obj.Unw.y(:,1:Ny),N,1), 0).*dMupyNw;
            vecY2 = max(reshape(obj.Uw.y(:,2:Ny+1),N,1), 0).*dMupyW + max(reshape(obj.Unw.y(:,2:Ny+1),N,1), 0).*dMupyNw;
            acc = pv/dt .* ( x1(:,1) .* obj.drho(:,1) .* s + x1(:,2) .* obj.drho(:,2) .* (1-s));
            DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
            DiagIndx = [-Nx, -1, 0, 1, Nx];
            J1p = J1p + spdiags(DiagVecs, DiagIndx, N, N);
            
            %% 2. Component 2 pressure block
            
            % 2.a divergence
            J2p = obj.T2;
            
            % 2.b: compressibility part
            dMupxW = obj.UpWindW.x*(obj.Mob(:, 1) .* x2(:,1) .* obj.drho(:,1));
            dMupyW = obj.UpWindW.y*(obj.Mob(:, 1) .* x2(:,1) .* obj.drho(:,1));
            dMupxNw = obj.UpWindNw.x*(obj.Mob(:, 2) .* x2(:,2) .* obj.drho(:,2));
            dMupyNw = obj.UpWindNw.y*(obj.Mob(:, 2) .* x2(:,2) .* obj.drho(:,2));
            
            vecX1 = min(reshape(obj.Uw.x(1:Nx,:),N,1), 0).*dMupxW + min(reshape(obj.Unw.x(1:Nx,:),N,1),0).*dMupxNw;
            vecX2 = max(reshape(obj.Uw.x(2:Nx+1,:),N,1), 0).*dMupxW + max(reshape(obj.Unw.x(2:Nx+1,:),N,1),0).*dMupxNw;
            vecY1 = min(reshape(obj.Uw.y(:,1:Ny),N,1), 0).*dMupyW + min(reshape(obj.Unw.y(:,1:Ny),N,1),0).*dMupyNw;
            vecY2 = max(reshape(obj.Uw.y(:,2:Ny+1),N,1), 0).*dMupyW + max(reshape(obj.Unw.y(:,2:Ny+1),N,1),0).*dMupyNw;
            acc = pv/dt .* (x2(:,1) .* obj.drho(:,1) .* s + x2(:,2) .* obj.drho(:,2) .* (1-s));
            DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
            DiagIndx = [-Nx, -1, 0, 1, Nx];
            J2p = J2p + spdiags(DiagVecs, DiagIndx, N, N);
            
            %% 3. Component 1 saturation block
            dMupxW = obj.UpWindW.x*(obj.dMob(:,1) .* x1(:,1) .* rho(:,1));
            dMupyW = obj.UpWindW.y*(obj.dMob(:,1) .* x1(:,1) .* rho(:,1));
            dMupxNw = obj.UpWindNw.x*(obj.dMob(:,2) .* x1(:,2) .* rho(:,2));
            dMupyNw = obj.UpWindNw.y*(obj.dMob(:,2) .* x1(:,2) .* rho(:,2));
            
            vecX1 = min(reshape(obj.Uw.x(1:Nx,:),N,1), 0).*dMupxW + min(reshape(obj.Unw.x(1:Nx,:),N,1), 0).*dMupxNw;
            vecX2 = max(reshape(obj.Uw.x(2:Nx+1,:),N,1), 0).*dMupxW + max(reshape(obj.Unw.x(2:Nx+1,:),N,1), 0).*dMupxNw;
            vecY1 = min(reshape(obj.Uw.y(:,1:Ny),N,1), 0).*dMupyW + min(reshape(obj.Unw.y(:,1:Ny),N,1), 0).*dMupyNw;
            vecY2 = max(reshape(obj.Uw.y(:,2:Ny+1),N,1), 0).*dMupyW + max(reshape(obj.Unw.y(:,2:Ny+1),N,1), 0).*dMupyNw;
            acc = pv/dt .* (x1(:,1) .* rho(:,1) - x1(:,2) .* rho(:,2));
            DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
            DiagIndx = [-Nx, -1, 0, 1, Nx];
            J1S = spdiags(DiagVecs,DiagIndx, N, N);
            
            %% 4. Component 2 saturation block
            dMupxW = obj.UpWindW.x*(obj.dMob(:, 1) .* x2(:,1) .* rho(:,1));
            dMupyW = obj.UpWindW.y*(obj.dMob(:, 1) .* x2(:,1) .* rho(:,1));
            dMupxNw = obj.UpWindNw.x*(obj.dMob(:, 2) .* x2(:,2) .* rho(:,2));
            dMupyNw = obj.UpWindNw.y*(obj.dMob(:, 2) .* x2(:,2) .* rho(:,2));
            
            vecX1 = min(reshape(obj.Uw.x(1:Nx,:),N,1), 0).*dMupxW + min(reshape(obj.Unw.x(1:Nx,:),N,1), 0).*dMupxNw;
            vecX2 = max(reshape(obj.Uw.x(2:Nx+1,:),N,1), 0).*dMupxW + max(reshape(obj.Unw.x(2:Nx+1,:),N,1), 0).*dMupxNw;
            vecY1 = min(reshape(obj.Uw.y(:,1:Ny),N,1), 0).*dMupyW + min(reshape(obj.Unw.y(:,1:Ny),N,1), 0).*dMupyNw;
            vecY2 = max(reshape(obj.Uw.y(:,2:Ny+1),N,1), 0).*dMupyW + max(reshape(obj.Unw.y(:,2:Ny+1),N,1), 0).*dMupyNw;
            acc = pv/dt .* (x2(:,1) .* rho(:,1) - x2(:,2) .* rho(:,2));
            DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
            DiagIndx = [-Nx, -1, 0, 1, Nx];
            J2S = spdiags(DiagVecs, DiagIndx, N, N);
            
            %% 5.Add capillarity
            %J1S = J1S - TMatrixW * spdiags(x1(:,1).* obj.dPc, 0, N, N);
            %J2S = J2S - TMatrixW * spdiags(x2(:,1).* obj.dPc, 0, N, N);
            
            %% Add wells to each block
            [J1p, J2p, J1S, J2S] = obj.AddWellsToJacobian(J1p, J2p, J1S, J2S, ProductionSystem.Reservoir.State, ProductionSystem.Wells, ProductionSystem.Reservoir.K);
            
            %% Full Jacobean: combine the 4 blocks
            Jacobian = [ J1p, J1S ;...
                         J2p, J2S ];
        end
        function delta = UpdateState(obj, delta, Status, FluidModel)
            % Update Solution
            Status.p = Status.p + delta(1:end/2);
            Status.S = Status.S + delta(end/2+1:end);
            
            % Remove unphysical solutions
            Status.S = min(Status.S,1);
            Status.S = max(Status.S,0);
            
            % Update density
            for i=1:FluidModel.NofPhases
                Status.rho(:, i) = FluidModel.Phases(i).ComputeDensity(Status.p);
            end
            
            % Update z
            Status.z = FluidModel.ComputeMassFractions(Status.S, Status.x1, Status.rho);
           
            %% Composition update loop
            Converged = 0;
            itCount = 1;
            s_0 = Status.S;
            while Converged == 0 && itCount < FluidModel.FlashSettings.MaxIt
                
                % 1. Stores old values
                s_old = Status.S;
                x_old = Status.x1;
                z_old = Status.z;
                
                % 2. Update Composition of the phases (Flash)
                SinglePhase = FluidModel.Flash(Status);
                
                % 3. Update x and S based on components mass balance
                FluidModel.ComputePhaseSaturation(Status, SinglePhase);
                
                % 4. Compute new total mass fractions (z)
                Status.z =  FluidModel.ComputeMassFractions(Status.S, Status.x1, Status.rho);
                
                % 5.a Compute errors
                InnerError1 = norm(abs(Status.S(:,1) - s_old(:,1)), inf);   %Checks if this loop is converging
                InnerError2 = norm(abs(Status.x1(:,2) - x_old(:,2)), inf);
                InnerError3 = norm(abs(Status.z(:,1) - z_old(:,1)), inf);
                
                % 5.b Check convergence
                if(InnerError1 < FluidModel.FlashSettings.TolInner && InnerError2 < FluidModel.FlashSettings.TolInner && InnerError3 < FluidModel.FlashSettings.TolInner)
                    Converged = 1;
                end
                
                itCount = itCount + 1;
            end
            delta_flash = Status.S - s_0;
            
            % Recompute delta
            delta(end/2+1:end) = delta(end/2+1:end) + delta_flash;
            
            % Update total density
            Status.rhoT = FluidModel.ComputeTotalDensity(Status.S, Status.rho);
            
            % Update Pc
            Status.pc = FluidModel.ComputePc(Status.S);
        end
        function  T = TransmissibilityMatrix(obj, Grid, Rho, x)
            %%%Transmissibility matrix construction
            Tx = zeros(Grid.Nx+1, Grid.Ny);
            Ty = zeros(Grid.Nx, Grid.Ny+1);
            
            %Apply upwind operator
            Mupx = obj.UpWindW.x*(obj.Mob(:,1) .* Rho(:,1) .* x(:,1)) + obj.UpWindNw.x*(obj.Mob(:,2) .* Rho(:,2) .* x(:,2));
            Mupy = obj.UpWindW.y*(obj.Mob(:,1) .* Rho(:,1) .* x(:,1)) + obj.UpWindNw.y*(obj.Mob(:,2) .* Rho(:,2) .* x(:,2));
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
        function [J1p, J2p, J1S, J2S] = AddWellsToJacobian(obj, J1p, J2p, J1S, J2S, State, Wells, K)
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
                end
                
            end
        end
        function G = ComputeGravityTerm(obj, N)
            G = zeros(N,N);
        end
    end
end