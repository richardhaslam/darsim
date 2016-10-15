% Overall Composition Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 15 October 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Overall_Composition_formulation < Compositional_formulation
    properties
        dxdz
    end
    methods
        function obj = Overall_Composition_formulation(n_components)
            obj@Compositional_formulation(n_components);
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State.S);
            obj.dMob = FluidModel.DMobDz(ProductionSystem.Reservoir.State);
            obj.drho = FluidModel.DrhoDp(ProductionSystem.Reservoir.State.p);
            obj.dPc = FluidModel.DPcDz(ProductionSystem.Reservoir.State);
            obj.dxdz = FluidModel.DxDz(ProductionSystem.Reservoir.State); % This is the bitchy part!! 
         end
        function Residual = BuildResidual(obj, DiscretizationModel, ProductionSystem)
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
            x(:,1:2) = ProductionSystem.Reservoir.State.x1;
            x(:,3:4) = 1 - x(:,1:2);
            
            % Density
            rho = ProductionSystem.Reservoir.State.rho;
            rhot = ProductionSystem.Reservoir.State.rhot;
            rhot_old = State0.rhot;
            
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
                n = rhot .* z(:,i);
                n_old = rhot_old .* z_old(:,i);
                % Phase Transmissibilities
                obj.TransmissibilityMatrix(DiscretizationModel.ReservoirGrid, rho, obj.GravityModel.RhoInt, x(:,(i-1)*2+1:(i-1)*2+2), i); 
                Residual = A * (n - n_old) ...             % Accumulation term
                           + obj.Tph{i, 1} *  P_ph1 ...    % Convective term                
                           + obj.Tph{i, 2} *  P_ph2...
                           + obj.Gph{i,1} * depth...       % Gravity
                           + obj.Gph{i,2} * depth...
                           - q(:,i);                       % Source/Sink  
            end
        end
        function Jacobian = BuildJacobian(obj)
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
            
            % Initialise local variables
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            N = DiscretizationModel.ReservoirGrid.N;
            %
            pv = DiscretizationModel.ReservoirGrid.Volume*ProductionSystem.Reservoir.Por;
            z = ProductionSystem.Reservoir.State.z;
            x(:,1:2) = ProductionSystem.Reservoir.State.x1;
            x(:,3:4) = 1 - x(:,1:2);
            rho = ProductionSystem.Reservoir.State.rho;
            rhot = ProductionSystem.Reservoir.State.rhot;
            
            % Fill in block by block
            Jp = cell(obj.NofComponents, 1);
            Jz = cell(obj.NofComponents, 1);
            for i=1:obj.NofComponents
                %% 1. Component i pressure block
                % 1.a: divergence
                Jp{i} = obj.Tph{i,1}  + obj.Tph{i, 2};
                % 1.b: compressibility part
                dMupxPh1 = obj.UpWind(1).x * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drho(:,1));
                dMupyPh1 = obj.UpWind(1).y * (obj.Mob(:, 1) .* x(:,(i-1)*2+1) .* obj.drho(:,1));
                dMupxPh2 = obj.UpWind(2).x * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drho(:,2));
                dMupyPh2 = obj.UpWind(2).y * (obj.Mob(:, 2) .* x(:,(i-1)*2+2) .* obj.drho(:,2));
                
                vecX1 = min(reshape(obj.U(1).x(1:Nx,:),N,1), 0).*dMupxPh1 + min(reshape(obj.U(2).x(1:Nx,:),N,1), 0).*dMupxPh2;
                vecX2 = max(reshape(obj.U(1).x(2:Nx+1,:),N,1), 0).*dMupxPh1 + max(reshape(obj.U(2).x(2:Nx+1,:),N,1), 0).*dMupxPh2;
                vecY1 = min(reshape(obj.U(1).y(:,1:Ny),N,1), 0).*dMupyPh1 + min(reshape(obj.U(2).y(:,1:Ny),N,1), 0).*dMupyPh2;
                vecY2 = max(reshape(obj.U(1).y(:,2:Ny+1),N,1), 0).*dMupyPh1 + max(reshape(obj.U(2).y(:,2:Ny+1),N,1), 0).*dMupyPh2;
                acc = pv/dt .* (obj.drhot .* z(:,i));
                DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
                DiagIndx = [-Nx, -1, 0, 1, Nx];
                Jp{i} = Jp{i} + spdiags(DiagVecs, DiagIndx, N, N);
                
                %% 2. Component i composition block
                dMupxPh1 = obj.UpWind(1).x*(obj.dMob(:,1) .* x(:,(i-1)*2+1) .* rho(:,1)...
                    + obj.Mob(:,1) .* obj.dxdz(:,(i-1)*2+1) .* rho(:,1));
                dMupyPh1 = obj.UpWind(1).y*(obj.dMob(:,1) .* x(:,(i-1)*2+1) .* rho(:,1)...
                    + obj.Mob(:,1) .* obj.dxdz(:,(i-1)*2+1) .* rho(:,1));
                dMupxPh2 = obj.UpWind(2).x*(obj.dMob(:,2) .* x(:,(i-1)*2+2) .* rho(:,2)...
                    + obj.Mob(:,2) .* obj.dxdz(:,(i-1)*2+2) .* rho(:,2));
                dMupyPh2 = obj.UpWind(2).y*(obj.dMob(:,2) .* x(:,(i-1)*2+2) .* rho(:,2)...
                    + obj.Mob(:,2) .* obj.dxdz(:,(i-1)*2+2) .* rho(:,2));
                
                vecX1 = min(reshape(obj.U(1).x(1:Nx,:),N,1), 0).*dMupxPh1 + min(reshape(obj.U(2).x(1:Nx,:),N,1), 0).*dMupxPh2;
                vecX2 = max(reshape(obj.U(1).x(2:Nx+1,:),N,1), 0).*dMupxPh1 + max(reshape(obj.U(2).x(2:Nx+1,:),N,1), 0).*dMupxPh2;
                vecY1 = min(reshape(obj.U(1).y(:,1:Ny),N,1), 0).*dMupyPh1 + min(reshape(obj.U(2).y(:,1:Ny),N,1), 0).*dMupyPh2;
                vecY2 = max(reshape(obj.U(1).y(:,2:Ny+1),N,1), 0).*dMupyPh1 + max(reshape(obj.U(2).y(:,2:Ny+1),N,1), 0).*dMupyPh2;
                acc = pv/dt .* rhot;
                DiagVecs = [-vecY2, -vecX2, vecY2+vecX2-vecY1-vecX1+acc, vecX1, vecY1];
                DiagIndx = [-Nx, -1, 0, 1, Nx];
                Jz{i} = spdiags(DiagVecs,DiagIndx, N, N);
                % Capillarity
                Jz{i} = Jz{i} - obj.Tph{i,1} * spdiags(obj.dPc, 0, N, N);               
            end
            
            %% 3. Add wells to each block
            [Jp{1}, Jp{2}, Jz{1}, Jz{2}] = ...
                obj.AddWellsToJacobian(Jp{1}, Jp{2}, Jz{1}, Jz{2},...
                ProductionSystem.Reservoir.State, ProductionSystem.Wells, ProductionSystem.Reservoir.K);
            
            %% Full Jacobian
            Jacobian = [ Jp{1}, Jz{1},...
                         Jp{2}, Jz{2}];
        end
        function UpdateState(obj, Status, delta, FluidModel)
            %% 1. Update Primary variables 
            Status.p = delta(1:end/2);
            Status.z(:,1) = delta(end/2+1:end);
            
            %% 2. Perform composition update
            SinglePhase = FluidModel.Flash(Status.z, Status.x1);
            
            %% 3. Compute Saturations
            FluidModel.ComputePhaseSaturation(Status.z, Status.x, SinglePhase);
            
            %% 4. Compute Densities
            FluidModel.ComputePhaseDensities(Status);
            
            %% 5. Compute Total Density
            Status.rhot = FluidModel.ComputeTotalDensity(obj, Status.S, Status.rho);
            
            %% 6. Compute Pc
            Status.pc = FluidModel.ComputePc(Status.S);   
        end
    end
end