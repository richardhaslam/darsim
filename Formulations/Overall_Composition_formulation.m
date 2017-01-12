% Overall Composition Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 16 December 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Overall_Composition_formulation < Compositional_formulation
    properties
        dxdz
        dxdp
        drhoTdp
        drhodz
        drhoTdz
        dMobdp
        dMobdz
    end
    methods
        function obj = Overall_Composition_formulation(n_components)
            obj@Compositional_formulation(n_components);
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.S);
            [obj.dxdp, obj.dxdz] = FluidModel.DxDpDz(ProductionSystem.Reservoir.State, obj.SinglePhase); % This is the bitchy part!! 
            obj.drhodp = FluidModel.DrhoDp(ProductionSystem.Reservoir.State.p);
            obj.drhodz = FluidModel.DrhoDz(ProductionSystem.Reservoir.State.z, obj.dxdz);
            dSdp = FluidModel.DSDp(ProductionSystem.Reservoir.State, obj.drhodp, -obj.dxdp(:,5));
            dSdz = FluidModel.DSDz(ProductionSystem.Reservoir.State, -obj.dxdz(:, end, :), obj.dxdz(:, 1));
            obj.drhoTdz = FluidModel.DrhotDz(ProductionSystem.Reservoir.State, obj.drhodz, dSdz);
            obj.drhoTdp = FluidModel.DrhotDp(ProductionSystem.Reservoir.State,obj.drhodp, dSdp);
            %obj.dMobdp = FluidModel.DMobDp(ProductionSystem.Reservoir.State, dSdp);
            obj.dMob = FluidModel.DMobDz(ProductionSystem.Reservoir.State, dSdz);
            obj.dPc = FluidModel.DPcDz(ProductionSystem.Reservoir.State, dSdz);
         end
        function Residual = BuildResidual(obj,ProductionSystem, DiscretizationModel, dt, State0)
            %Create local variables
            N = DiscretizationModel.ReservoirGrid.N;
          
            % Pore volume
            pv = ProductionSystem.Reservoir.Por*DiscretizationModel.ReservoirGrid.Volume;
            
            % Phase pressures
            P_ph1 = ProductionSystem.Reservoir.State.p - ProductionSystem.Reservoir.State.pc;
            P_ph2 = ProductionSystem.Reservoir.State.p;
            % overall compositions
            z_old = State0.z;
            z = ProductionSystem.Reservoir.State.z;
            
            % Compositions 
            x = ProductionSystem.Reservoir.State.x;
            
            % Density
            rho = ProductionSystem.Reservoir.State.rho;
            rhoT = ProductionSystem.Reservoir.State.rhoT;
            rhoT_old = State0.rhoT;
            
            % Depths
            depth = DiscretizationModel.ReservoirGrid.Depth;
            
            %Accumulation term
            A = speye(N)*pv/dt;
            
            %Source terms
            q = obj.ComputeSourceTerms(N, ProductionSystem.Wells);
            
            %%  Build Residual
            Residual = zeros(N*obj.NofComponents, 1);
            for i=1:obj.NofComponents
                % Total moles
                n = rhoT .* z(:,i);
                n_old = rhoT_old .* z_old(:,i);
                % Phase Transmissibilities
                obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, rho, obj.GravityModel.RhoInt, x(:,(i-1)*2+1:(i-1)*2+2), i); 
                Residual((i-1)*N+1:i*N) = A * (n - n_old) ...             % Accumulation term
                           + obj.Tph{i, 1} *  P_ph1 ...    % Convective term                
                           + obj.Tph{i, 2} *  P_ph2...
                           + obj.Gph{i,1} * depth...       % Gravity
                           + obj.Gph{i,2} * depth...
                           - q(:,i);                       % Source/Sink  
            end
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            
            % Initialise local variables
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N = DiscretizationModel.ReservoirGrid.N;
            %
            pv = DiscretizationModel.ReservoirGrid.Volume*ProductionSystem.Reservoir.Por;
            z = ProductionSystem.Reservoir.State.z;
            x = ProductionSystem.Reservoir.State.x;
            rho = ProductionSystem.Reservoir.State.rho;
            rhoT = ProductionSystem.Reservoir.State.rhoT;
            
            % Fill in block by block
            Jp = cell(obj.NofComponents, 1);
            Jz = cell(obj.NofComponents, obj.NofComponents-1);
            for i=1:obj.NofComponents
                %% 1. Component i pressure block
                % 1.a: divergence
                Jp{i} = obj.Tph{i,1}  + obj.Tph{i, 2};
                % 1.b: compressibility part
                dMupxPh1 = obj.UpWind(1).x * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drhodp(:,1) + obj.Mob(:, 1) .* obj.dxdp(:,(i-1)*2+1) .* rho(:,1));
                dMupyPh1 = obj.UpWind(1).y * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drhodp(:,1) + obj.Mob(:, 1) .* obj.dxdp(:,(i-1)*2+1) .* rho(:,1));
                dMupzPh1 = obj.UpWind(1).z * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drhodp(:,1) + obj.Mob(:, 1) .* obj.dxdp(:,(i-1)*2+1) .* rho(:,1));
                dMupxPh2 = obj.UpWind(2).x * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drhodp(:,2) + obj.Mob(:, 2) .* obj.dxdp(:,(i-1)*2+2) .* rho(:,2));
                dMupyPh2 = obj.UpWind(2).y * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drhodp(:,2) + obj.Mob(:, 2) .* obj.dxdp(:,(i-1)*2+2) .* rho(:,2));
                dMupzPh2 = obj.UpWind(2).z * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drhodp(:,2) + obj.Mob(:, 2) .* obj.dxdp(:,(i-1)*2+2) .* rho(:,2));
                
                vecX1 = min(reshape(obj.U(1).x(1:Nx,:,:),N,1), 0).*dMupxPh1 + min(reshape(obj.U(2).x(1:Nx,:,:),N,1), 0).*dMupxPh2;
                vecX2 = max(reshape(obj.U(1).x(2:Nx+1,:,:),N,1), 0).*dMupxPh1 + max(reshape(obj.U(2).x(2:Nx+1,:,:),N,1), 0).*dMupxPh2;
                vecY1 = min(reshape(obj.U(1).y(:,1:Ny,:),N,1), 0).*dMupyPh1 + min(reshape(obj.U(2).y(:,1:Ny,:),N,1), 0).*dMupyPh2;
                vecY2 = max(reshape(obj.U(1).y(:,2:Ny+1,:),N,1), 0).*dMupyPh1 + max(reshape(obj.U(2).y(:,2:Ny+1,:),N,1), 0).*dMupyPh2;
                vecZ1 = min(reshape(obj.U(1).z(:,:,1:Nz),N,1), 0).*dMupzPh1 + min(reshape(obj.U(2).z(:,:,1:Nz),N,1), 0).*dMupzPh2;
                vecZ2 = max(reshape(obj.U(1).z(:,:,2:Nz+1),N,1), 0).*dMupzPh1 + max(reshape(obj.U(2).z(:,:,2:Nz+1),N,1), 0).*dMupzPh2;
                acc = pv/dt .* (obj.drhoTdp .* z(:,i));
                DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2 + vecY2+vecX2-vecY1-vecX1-vecZ1+acc, vecX1, vecY1, vecZ1];
                DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                Jp{i} = Jp{i} + spdiags(DiagVecs, DiagIndx, N, N);
                
                %% 2. Component i composition block
                for j=1:obj.NofComponents-1
                    dMupxPh1 = obj.UpWind(1).x*...
                        ( obj.dMob(:, 1, j) .* x(:,(i-1)*2+1) .* rho(:,1)...
                        + obj.Mob(:, 1) .* obj.dxdz(:,(i-1)*2+1, j) .* rho(:,1)...
                        + obj.Mob(:, 1) .* x(:, (i-1)*2+1) .* obj.drhodz(:, 1));
                    dMupyPh1 = obj.UpWind(1).y*...
                        ( obj.dMob(:, 1, j) .* x(:, (i-1)*2+1) .* rho(:,1)...
                        + obj.Mob(:, 1) .* obj.dxdz(:, (i-1)*2+1, j) .* rho(:,1)...
                        + obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drhodz(:, 1));
                    dMupzPh1 = obj.UpWind(1).z*...
                        ( obj.dMob(:, 1, j) .* x(:,(i-1)*2+1) .* rho(:,1)...
                        + obj.Mob(:,1) .* obj.dxdz(:,(i-1)*2+1) .* rho(:,1)...
                        + obj.Mob(:,1) .* x(:,(i-1)*2+1) .* obj.drhodz(:, 1));
                    dMupxPh2 = obj.UpWind(2).x*...
                        ( obj.dMob(:, 2, j) .* x(:,(i-1)*2+2) .* rho(:,2)...
                        + obj.Mob(:,2) .* obj.dxdz(:,(i-1)*2+2) .* rho(:,2)...
                        + obj.Mob(:,2) .* x(:,(i-1)*2+2) .* obj.drhodz(:, 2));
                    dMupyPh2 = obj.UpWind(2).y*...
                        ( obj.dMob(:, 2, j) .* x(:,(i-1)*2+2) .* rho(:,2)...
                        + obj.Mob(:, 2) .* obj.dxdz(:,(i-1)*2+2, j) .* rho(:,2)...
                        + obj.Mob(:, 2) .* x(:, (i-1)*2+2) .* obj.drhodz(:, 2));
                    dMupzPh2 = obj.UpWind(2).z*...
                        ( obj.dMob(:, 2, j) .* x(:, (i-1)*2+2) .* rho(:,2)...
                        + obj.Mob(:,2) .* obj.dxdz(:, (i-1)*2+2, j) .* rho(:,2)...
                        + obj.Mob(:,2) .* x(:, (i-1)*2+2) .* obj.drhodz(:, 2));
                    
                    vecX1 = min(reshape(obj.U(1).x(1:Nx,:,:),N,1), 0).*dMupxPh1 + min(reshape(obj.U(2).x(1:Nx,:,:),N,1), 0).*dMupxPh2;
                    vecX2 = max(reshape(obj.U(1).x(2:Nx+1,:,:),N,1), 0).*dMupxPh1 + max(reshape(obj.U(2).x(2:Nx+1,:,:),N,1), 0).*dMupxPh2;
                    vecY1 = min(reshape(obj.U(1).y(:,1:Ny,:),N,1), 0).*dMupyPh1 + min(reshape(obj.U(2).y(:,1:Ny,:),N,1), 0).*dMupyPh2;
                    vecY2 = max(reshape(obj.U(1).y(:,2:Ny+1,:),N,1), 0).*dMupyPh1 + max(reshape(obj.U(2).y(:,2:Ny+1,:),N,1), 0).*dMupyPh2;
                    vecZ1 = min(reshape(obj.U(1).z(:,:,1:Nz),N,1), 0).*dMupzPh1 + min(reshape(obj.U(2).z(:,:,1:Nz),N,1), 0).*dMupzPh2;
                    vecZ2 = max(reshape(obj.U(1).z(:,:,2:Nz+1),N,1), 0).*dMupzPh1 + max(reshape(obj.U(2).z(:,:,2:Nz+1),N,1), 0).*dMupzPh2;
                    
                    if i == obj.NofComponents
                        acc = - pv/dt .* rhoT + pv/dt .* z(:,i) .* obj.drhoTdz(:, j);
                    else
                        acc = pv/dt .* rhoT + pv/dt .* z(:,i) .* obj.drhoTdz(:, j);
                    end
                    
                    DiagVecs = [-vecZ2, -vecY2, -vecX2, vecZ2+vecY2+vecX2-vecZ1-vecY1-vecX1+acc, vecX1, vecY1, vecZ1];
                    DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
                    Jz{i,j} = spdiags(DiagVecs,DiagIndx, N, N);
                    % Capillarity
                    Jz{i,j} = Jz{i,j} - obj.Tph{i,1} * spdiags(obj.dPc, 0, N, N);
                end
            end
            
            %% 3. Add wells to each block
            [Jp, Jz] = obj.AddWellsToJacobian(Jp, Jz, ProductionSystem.Reservoir.State, ProductionSystem.Wells, ProductionSystem.Reservoir.K);
            
            %% Full Jacobian
            if obj.NofComponents <= 2 
            Jacobian = [Jp{1}, Jz{1,1}
                        Jp{2}, Jz{2,1}];
            else
            % for now let's consider a max of 3 components!! 
            Jacobian = [Jp{1}, Jz{1,1}, Jz{1,2};...
                        Jp{2}, Jz{2,1}, Jz{2,2};...
                        Jp{3}, Jz{3,1}, Jz{3,2}];
            end
        end
        function [Jp, Jz] = AddWellsToJacobian(obj, Jp, Jz, State, Wells, K)
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %Injectors
            for i=1:Wells.NofInj
                a = Inj(i).Cells;
                for j=1:length(a)
                    Jp{1}(a(j),a(j)) = Jp{1}(a(j),a(j)) ...
                                       + Inj(i).PI * K(a(j)) * (Inj(i).Mob(:,1) * Inj(i).rho(j,1) * Inj(i).x(:,1) ...
                                       + Inj(i).Mob(:,2) *  Inj(i).rho(j,2) * Inj(i).x(:,2));
                    Jp{2}(a(j),a(j)) = Jp{2}(a(j),a(j)) ...
                                       + Inj(i).PI * K(a(j)) * (Inj(i).Mob(:,1) * Inj(i).rho(j,1) * Inj(i).x(:,3) ...
                                       + Inj(i).Mob(:,2) *  Inj(i).rho(j,2) * Inj(i).x(:,4));
                end
            end
            
            % Compositions 
            x = State.x;
            % Density
            rho = State.rho;
            
            %Producers
            for i=1:Wells.NofProd
                b = Prod(i).Cells;
                for j=1:length(b)
                    %Pressure blocks
                    Jp{1}(b(j),b(j)) = Jp{1}(b(j),b(j))...
                        + Prod(i).PI * K(b(j)) *...
                        ( obj.Mob(b(j), 1) * State.rho(b(j),1) * x(b(j),1) ...
                        + obj.Mob(b(j), 2) * rho(b(j),2) * x(b(j),2)...
                        )...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * x(b(j),1) * obj.drhodp(b(j),1) * (Prod(i).p(j) - State.p(b(j))) ...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * obj.dxdp(b(j), 1) * rho(b(j),1) * (Prod(i).p(j) - State.p(b(j))) ...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * x(b(j),2) * obj.drhodp(b(j),2) * (Prod(i).p(j) - State.p(b(j)))...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * obj.dxdp(b(j),2) * rho(b(j),2) * (Prod(i).p(j) - State.p(b(j)));
                    Jp{2}(b(j),b(j)) = Jp{2}(b(j),b(j)) ...
                        + Prod(i).PI * K(b(j)) * (obj.Mob(b(j), 1) * rho(b(j),1) * x(b(j), 3) ...
                        + obj.Mob(b(j), 2) * rho(b(j),2) * x(b(j), 4))...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * x(b(j),3) * obj.drhodp(b(j), 1) * (Prod(i).p(j) - State.p(b(j))) ...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * obj.dxdp(b(j),3) * rho(b(j), 1) * (Prod(i).p(j) - State.p(b(j))) ...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * x(b(j),3) * obj.drhodp(b(j), 2) * (Prod(i).p(j) - State.p(b(j))) ...
                        - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * obj.dxdp(b(j),4) * rho(b(j), 2) * (Prod(i).p(j) - State.p(b(j)));
                    % Overall composition blocks
                    Jz{1}(b(j),b(j)) = Jz{1}(b(j),b(j))...
                        - Prod(i).PI * K(b(j)) *...
                        ( obj.dMob(b(j), 1) * rho(b(j), 1) * x(b(j), 1) ...
                        + obj.dMob(b(j), 2) * rho(b(j), 2) * x(b(j), 2)...
                        + obj.Mob(b(j), 1) * rho(b(j), 1) * obj.dxdz(b(j), 1) ...
                        + obj.Mob(b(j), 2) * rho(b(j), 2) * obj.dxdz(b(j), 2)...
                        ) * (Prod(i).p(j) - State.p(b(j)));
                    Jz{2}(b(j),b(j)) = Jz{2}(b(j),b(j)) ...
                        - Prod(i).PI * K(b(j)) *...
                        ( obj.dMob(b(j), 1) * State.rho(b(j),1) * x(b(j), 3) ...
                        + obj.dMob(b(j), 2) * rho(b(j),2) * x(b(j), 4)...
                        + obj.Mob(b(j), 1) * rho(b(j),1) * obj.dxdz(b(j), 3) ...
                        + obj.Mob(b(j), 2) * rho(b(j),2) * obj.dxdz(b(j), 4)...
                        ) * (Prod(i).p(j) - State.p(b(j)));
                end
                
            end
        end
        function UpdateState(obj, delta, Status, FluidModel)
            %% 1. Update Primary variables 
            Status.p = Status.p + delta(1:end/2);
            Status.z(:,1) = Status.z(:,1) + delta(end/2+1:end);
            Status.z(:,2) = 1 - Status.z(:,1);           
            % Update remaining variables 
            if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                return
            else
                obj.UpdatePhaseCompositions(Status, FluidModel);     
            end
        end
    end
end