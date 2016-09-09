% Immiscible Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Immiscible_formulation < fim_formulation
    properties
        Tw
        Tnw
        Mob
        dMob
        dPc
    end
    methods
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.S);
            obj.dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State.S);
            obj.dPc = FluidModel.DPcDS(ProductionSystem.Reservoir.State.S);
        end
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            %Create local variables
            N = DiscretizationModel.ReservoirGrid.N;
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            p_old = State0.p;
            s_old = State0.S;
            rho_old = State0.rho;
            p = ProductionSystem.Reservoir.State.p;
            s = ProductionSystem.Reservoir.State.S;
            rho = ProductionSystem.Reservoir.State.rho;
            pc = ProductionSystem.Reservoir.State.pc;
            
            %Accumulation Term
            Ap = speye(N)*0; %It's still incompressible
            AS = speye(N)*pv/dt;
            
            %Transmissibility matrix
            obj.Tw = obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, obj.UpWindW, obj.Mob(:,1), rho(:, 1));
            obj.Tnw = obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, obj.UpWindNw, obj.Mob(:,2), rho(:, 2));
            
            %Gravity
            G = obj.ComputeGravityTerm(N);
            
            %Source terms
            [qw, qnw] = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            
            %% RESIDUAL
            %Wetting phase (This one has capillarity)
            Rw = Ap*(p-p_old) + AS*(rho(:,1) .* s - rho_old(:,1) .* s_old) + obj.Tw*p - obj.Tw*pc + G*s - qw;
            %Non-wetting phase
            Rnw = Ap*(p-p_old) - AS*(rho(:,2) .* s - rho_old(:,2) .* s_old) + obj.Tnw*p + G*s - qnw;
            
            %Stick them together
            Residual = [Rw; Rnw];
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            %Build FIM Jacobian
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            N = DiscretizationModel.ReservoirGrid.N;
            pv = DiscretizationModel.ReservoirGrid.Volume*ProductionSystem.Reservoir.Por;
            rho = ProductionSystem.Reservoir.State.rho;
            
            % BUILD FIM JACOBIAN BLOCK BY BLOCK         
            %1. Rnw Pressure Block
            Jnwp = obj.Tnw;
            
            %2. Rw Pressure Block
            Jwp = obj.Tw;
            
            %3. Rnw Saturation Block
            dMupxNw = obj.UpWindNw.x * (obj.dMob(:,2) .* rho(:, 2));
            dMupyNw = obj.UpWindNw.y * (obj.dMob(:,2) .* rho(:,2));
            %Construct Jnws block
            x1 = min(reshape(obj.Unw.x(1:Nx,:),N,1),0).*dMupxNw;
            x2 = max(reshape(obj.Unw.x(2:Nx+1,:),N,1),0).*dMupxNw;
            y1 = min(reshape(obj.Unw.y(:,1:Ny),N,1),0).*dMupyNw;
            y2 = max(reshape(obj.Unw.y(:,2:Ny+1),N,1),0).*dMupyNw;
            v = ones(N,1)*pv/dt .* rho(:,2);
            DiagVecs = [-y2, -x2, y2+x2-y1-x1-v, x1, y1];
            DiagIndx = [-Nx, -1, 0, 1, Nx];
            JnwS = spdiags(DiagVecs,DiagIndx,N,N);
            
            %4. Rw Saturation Block
            dMupxw = obj.UpWindW.x * (obj.dMob(:,1) .* rho(:,1));
            dMupyw = obj.UpWindW.y * (obj.dMob(:,1) .* rho(:,1));
            %Construct JwS block
            x1 = min(reshape(obj.Uw.x(1:Nx,:),N,1),0).*dMupxw;
            x2 = max(reshape(obj.Uw.x(2:Nx+1,:),N,1),0).*dMupxw;
            y1 = min(reshape(obj.Uw.y(:,1:Ny),N,1),0).*dMupyw;
            y2 = max(reshape(obj.Uw.y(:,2:Ny+1),N,1),0).*dMupyw;
            v = ones(N,1)*pv/dt .* rho(:,1);
            DiagVecs = [-y2, -x2, y2+x2-y1-x1+v, x1, y1];
            DiagIndx = [-Nx, -1, 0, 1, Nx];
            JwS = spdiags(DiagVecs,DiagIndx,N,N);
            
            %Add capillarity
            CapJwS = Jwp * spdiags(obj.dPc, 0, N, N);
            JwS = JwS - CapJwS;
            
            % Add Wells
            [Jwp, JwS, Jnwp, JnwS] = obj.AddWellsToJacobian(Jwp, JwS, Jnwp, JnwS, ProductionSystem.Reservoir.State, ProductionSystem.Wells, ProductionSystem.Reservoir.K(:,1));
            
            % Full Jacobian: put the 4 blocks together
            Jacobian = [Jwp, JwS; Jnwp, JnwS];
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
            
            % Update total density
            Status.rhoT = FluidModel.ComputeTotalDensity(Status.S, Status.rho);  
            
            % Update Pc
            Status.pc = FluidModel.ComputePc(Status.S);
        end
        function T = TransmissibilityMatrix(obj, Grid, UpWind, M, rho)
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            N = Grid.N;
            
            %%%Transmissibility matrix construction
            Tx = zeros(Nx+1, Ny);
            Ty = zeros(Nx, Ny+1);
            %Apply upwind operator
            Mupx = UpWind.x*(M.*rho);
            Mupy = UpWind.y*(M.*rho);
            Mupx = reshape(Mupx, Nx, Ny);
            Mupy = reshape(Mupy, Nx, Ny);
            Tx(2:Nx,:)= Grid.Tx(2:Nx,:).*Mupx(1:Nx-1,:);
            Ty(:,2:Ny)= Grid.Ty(:,2:Ny).*Mupy(:,1:Ny-1);
            %Construct matrix
            x1 = reshape(Tx(1:Nx,:),N,1);
            x2 = reshape(Tx(2:Nx+1,:),N,1);
            y1 = reshape(Ty(:,1:Ny),N,1);
            y2 = reshape(Ty(:,2:Ny+1),N,1);
            DiagVecs = [-y2,-x2,y2+x2+y1+x1,-x1,-y1];
            DiagIndx = [-Nx,-1,0,1,Nx];
            T = spdiags(DiagVecs,DiagIndx,N,N);
        end
        function [qw, qnw] = ComputeSourceTerms(obj, N, Wells)
            q = zeros(N, 2);    
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
            qw = q(:,1);
            qnw = q(:,2);
        end
        function [Jwp, JwS, Jnwp, JnwS] = AddWellsToJacobian(obj, Jwp, JwS, Jnwp, JnwS, State, Wells, K)
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                for j=1:length(a)
                    Jnwp(a(j),a(j)) = Jnwp(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mob(:, 2)*Inj.rho(:, 2);
                    Jwp(a(j),a(j)) = Jwp(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mob(:, 1)*Inj.rho(:, 1);
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                for j=1:length(b)
                    Jnwp(b(j),b(j)) = Jnwp(b(j),b(j)) + Prod(i).PI*K(b(j)).*obj.Mob(b(j), 2) .* State.rho(b(j), 2);
                    Jwp(b(j),b(j)) = Jwp(b(j),b(j)) + Prod(i).PI*K(b(j)).*obj.Mob(b(j), 1) .* State.rho(b(j), 1);
                    JwS(b(j),b(j)) = JwS(b(j),b(j)) - Prod(i).PI*K(b(j)).* State.rho(b(j), 1) .* (Prod(i).p - State.p(b(j))).*obj.dMob(b(j), 1);
                    JnwS(b(j),b(j)) = JnwS(b(j),b(j)) - Prod(i).PI*K(b(j)).* State.rho(b(j), 2) .* (Prod(i).p - State.p(b(j))).*obj.dMob(b(j), 2);
                end
            end
        end
        function G = ComputeGravityTerm(obj, N)
            G = zeros(N);
        end
    end
end