% Immiscible Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 26 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Immiscible_formulation < fim_formulation
    properties
    end
    methods
        function obj = Immiscible_formulation()
            obj@fim_formulation();
            obj.Tph = cell(2,1);
            obj.Gph = cell(2,1);
        end
        function Reset(obj)
        end
        function SavePhaseState(obj)
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            obj.drhodp = FluidModel.DrhoDp(ProductionSystem.Reservoir.State.p);
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.S);
            obj.dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State.S);
            obj.dPc = FluidModel.DPcDS(ProductionSystem.Reservoir.State.S);
        end
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Create local variables
            N = DiscretizationModel.ReservoirGrid.N;
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            s_old(:,1) = State0.S;
            s_old(:,2) = 1 - s_old(:,1);
            rho_old = State0.rho;
            s(:,1) = ProductionSystem.Reservoir.State.S;
            s(:,2) = 1 - s(:,1);
            rho = ProductionSystem.Reservoir.State.rho;
            P(:,1) = ProductionSystem.Reservoir.State.p - ProductionSystem.Reservoir.State.pc;
            P(:,2) = ProductionSystem.Reservoir.State.p;
            depth = DiscretizationModel.ReservoirGrid.Depth;
            
            % Accumulation Term
            AS = speye(N)*pv/dt;
            
            % Transmissibility matrix
            obj.TransmissibilityMatrix (DiscretizationModel.ReservoirGrid, rho, obj.GravityModel.RhoInt);
                  
            % Source terms
            [q(:,1), q(:,2)] = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            
            %% RESIDUAL
            Residual = zeros(2*N,1);
            for i=1:2
                Residual((i-1)*N+1:i*N)  = AS*(rho(:,i) .* s(:,i) - rho_old(:,i) .* s_old(:,i))...
                                           + obj.Tph{i} * P(:, i)...
                                           + obj.Gph{i} * depth...
                                           - q(:,i);
            end   
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            %Build FIM Jacobian
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N = DiscretizationModel.ReservoirGrid.N;
            pv = DiscretizationModel.ReservoirGrid.Volume*ProductionSystem.Reservoir.Por;
            rho = ProductionSystem.Reservoir.State.rho;
            s(:,1) = ProductionSystem.Reservoir.State.S;
            s(:,2) = 1 - s(:, 1);
            
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
                
                vecX1 = min(reshape(obj.U(i).x(1:Nx,:,:), N, 1), 0).*dMupx;
                vecX2 = max(reshape(obj.U(i).x(2:Nx+1,:,:), N, 1), 0).*dMupx;
                vecY1 = min(reshape(obj.U(i).y(:,1:Ny,:), N, 1), 0).*dMupy;
                vecY2 = max(reshape(obj.U(i).y(:,2:Ny+1,:), N, 1), 0).*dMupy;
                vecZ1 = min(reshape(obj.U(i).z(:,:,1:Nz), N, 1), 0).*dMupz;
                vecZ2 = max(reshape(obj.U(i).z(:,:,2:Nz+1), N, 1), 0).*dMupz; 
                acc = pv/dt .* obj.drhodp(:,i) .* s(:,i);
                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                Jp{i} = Jp{i} + spdiags(DiagVecs, DiagIndx, N, N);
                
                % 2. Saturation Block
                dMupx = obj.UpWind(i).x * (obj.dMob(:,i) .* rho(:,i));
                dMupy = obj.UpWind(i).y * (obj.dMob(:,i) .* rho(:,i));
                dMupz = obj.UpWind(i).z * (obj.dMob(:,i) .* rho(:,i));
                % Construct JS block
                x1 = min(reshape(obj.U(i).x(1:Nx,:,:), N, 1),0).*dMupx;
                x2 = max(reshape(obj.U(i).x(2:Nx+1,:,:), N, 1),0).*dMupx;
                y1 = min(reshape(obj.U(i).y(:,1:Ny,:), N, 1),0).*dMupy;
                y2 = max(reshape(obj.U(i).y(:,2:Ny+1,:), N, 1),0).*dMupy;
                z1 = min(reshape(obj.U(i).z(:,:,1:Nz), N, 1),0).*dMupz;
                z2 = max(reshape(obj.U(i).z(:,:,2:Nz+1), N, 1),0).*dMupz;
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
        function delta = UpdateState(obj, delta, Status, FluidModel)
            % Update Solution
            Status.p = Status.p + delta(1:end/2);
            Status.S = Status.S + delta(end/2+1:end);
            
            % Remove unphysical solutions
            Status.S = min(Status.S,1);
            Status.S = max(Status.S,0);
            
            % Update density
            FluidModel.ComputePhaseDensities(Status);
            
            % Update z
            Status.z = FluidModel.ComputeTotalFractions(Status.S, Status.x, Status.rho);
            
            % Update total density
            Status.rhoT = FluidModel.ComputeTotalDensity(Status.S, Status.rho);  
            
            % Update Pc
            Status.pc = FluidModel.ComputePc(Status.S);
        end
        function TransmissibilityMatrix(obj, Grid, rho, RhoInt)
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            
            for i=1:2
                % Transmissibility matrix construction
                Tx = zeros(Nx+1, Ny, Nz);
                Ty = zeros(Nx, Ny+1, Nz);
                Tz = zeros(Nx, Ny, Nz+1);
                
                % Apply upwind operator
                Mupx = obj.UpWind(i).x*(obj.Mob(:,i) .* rho(:,i));
                Mupy = obj.UpWind(i).y*(obj.Mob(:,i) .* rho(:,i));
                Mupz = obj.UpWind(i).z*(obj.Mob(:,i) .* rho(:,i));
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
                obj.Tph{i} = spdiags(DiagVecs,DiagIndx,N,N);
                
                % Gravity Matrix
                Tx(2:Grid.Nx,:,:)= Tx(2:Grid.Nx,:,:) .* RhoInt(i).x(2:Grid.Nx,:,:);
                Ty(:,2:Grid.Ny,:)= Ty(:,2:Grid.Ny,:) .* RhoInt(i).y(:,2:Grid.Ny,:);
                Tz(:,:,2:Grid.Nz)= Tz(:,:,2:Grid.Nz) .* RhoInt(i).z(:,:,2:Grid.Nz);
                
               % Construct matrix
                x1 = reshape(Tx(1:Nx,:,:), N, 1);
                x2 = reshape(Tx(2:Nx+1,:,:), N, 1);
                y1 = reshape(Ty(:,1:Ny,:), N, 1);
                y2 = reshape(Ty(:,2:Ny+1,:), N, 1);
                z1 = reshape(Tz(:,:,1:Nz), N, 1);
                z2 = reshape(Tz(:,:,2:Nz+1), N, 1);
                DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
                DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
                obj.Gph{i} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);     
            end    
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
                    Jnwp(a(j),a(j)) = Jnwp(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mob(:, 2)*Inj(i).rho(:, 2);
                    Jwp(a(j),a(j)) = Jwp(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mob(:, 1)*Inj(i).rho(:, 1);
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                for j=1:length(b)
                    Jnwp(b(j),b(j)) = Jnwp(b(j),b(j)) + Prod(i).PI*K(b(j)).*obj.Mob(b(j), 2) .* State.rho(b(j), 2)...
                     - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * obj.drhodp(b(j),2) * (Prod(i).p - State.p(b(j)));                    
                    Jwp(b(j),b(j)) = Jwp(b(j),b(j)) + Prod(i).PI*K(b(j)).*obj.Mob(b(j), 1) .* State.rho(b(j), 1)...
                     - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * obj.drhodp(b(j),1) * (Prod(i).p - State.p(b(j)));
                    
                    JwS(b(j),b(j)) = JwS(b(j),b(j)) - Prod(i).PI*K(b(j)).* State.rho(b(j), 1) .* (Prod(i).p - State.p(b(j))).*obj.dMob(b(j), 1);
                    JnwS(b(j),b(j)) = JnwS(b(j),b(j)) - Prod(i).PI*K(b(j)).* State.rho(b(j), 2) .* (Prod(i).p - State.p(b(j))).*obj.dMob(b(j), 2);
                end
            end
        end
    end
end