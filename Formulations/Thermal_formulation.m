% Thermal formulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadytio
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Thermal_formulation < formulation
    properties
        drho % make them 2 columns --> first column is drhodp secon column is drhodT
        dmu
        dh
    end
    methods
        function obj = Thermal_formulation()
            obj@formulation();
            obj.Tph = cell(2,1);
            obj.Gph = cell(2,1);
        end
        function x = GetPrimaryUnknowns(obj, ProductionSystem, DiscretizationModel)
             Nt = DiscretizationModel.N;
             Nm = DiscretizationModel.ReservoirGrid.N;
             x = zeros(obj.NofPhases * Nt, 1);
             %% Reservoir
             Start = 1;
             End = Nm;
             x(Start:End) = ProductionSystem.Reservoir.State.Properties('P_2').Value;
        end
        function x = GetPrimaryPressure(obj, ProductionSystem, DiscretizationModel)
            Nt = DiscretizationModel.N;
            Nm = DiscretizationModel.ReservoirGrid.N;
            %Nf = DiscretizationModel.FracturesGrid.N;
            x = zeros(Nt, 1);
            if obj.NofPhases > 1
                x(1:Nm) = ProductionSystem.Reservoir.State.Properties('P_2').Value;
            else
                x(1:Nm) = ProductionSystem.Reservoir.State.Properties('P_1').Value;
            end
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            %% 1. Reservoir Properteis and Derivatives
            obj.drho(:, 1) = FluidModel.DrhoDp(ProductionSystem.Reservoir.State);
            obj.drho(:, 2) = FluidModel.DrhoDT(ProductionSystem.Reservoir.State);
            % same thing for the other variables (mu and h)
            obj.Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State);
            obj.dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State);
        end
        %% Methods for FIM Coupling
        function Residual = BuilResidual(obj)
          % Here you build residual
        end
        function J = BuildJacobian(obj)
            % Jacobian construction
        end
        function delta = UpdateState(obj, delta, ProductionSystem, FluidModel, DiscretizationModel)
            if sum(isnan(delta))
                % if the solution makes no sense, skip this step
                return
            else
                N = DiscretizationModel.N;
                deltaP = delta(1:N);
                deltaT = delta(N+1: end);
                %% Update reservoir state
                % 1. Update Pressure
                Pm = ProductionSystem.Reservoir.State.Properties(['P_1']);
                Pm.update(deltaP);
                % 2. Update Temperature
                Tm = ProductionSystem.Reservoir.State.Properties('T');
                Tm.update(deltaT);
                
                % Update remaining properteis
                % 3. Update Phase Densities
                FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
                % 4. Update viscosity
                FluidModel.ComputePhaseViscosities(ProductionSystem.Reservoir.State);
                % 5. Update enthalpy
                FluidModel.ComputePhaseEnthalpies(ProductionSystem.Reservoir.State);
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
        function [Jp, JS] = AddWellsToJacobian(obj, Jp, JS, State, Wells, K, ph)
            % Define Local handles
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %p = State.Properties('P_2').Value;
            %rho = State.Properties(['rho_', num2str(ph)]).Value;
          
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                [dQdp, ~] = Inj(i).dQPhasesdPdS(K, obj.NofPhases);
                for j=1:length(a)
                    Jp(a(j),a(j)) = Jp(a(j),a(j)) - dQdp(j, ph);
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                [dQdp, dQdS] = Prod(i).dQPhasesdPdS(State, K, obj.Mob, obj.dMob, obj.drhodp, obj.NofPhases);
                for j=1:length(b)
                    Jp(b(j),b(j)) = Jp(b(j),b(j)) - dQdp(j, ph);                    
                    JS(b(j),b(j)) = JS(b(j),b(j)) - dQdS(j, ph);
                end
            end
        end
        function CFL = ComputeCFLNumber(obj, ProductionSystem, DiscretizationModel, dt)
            % Not important for now
            CFL = 0;
        end
    end
end