%  Sequential Formulation base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 22 July 2016
%Last modified: 26 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef seq_formulation < handle
    properties
        NumOfEquations
        Mob
        Mobt
        f
        df
        U
        Tx
        Ty
        Tz
        Qwells
        V
        GravityModel
        SinglePhase
    end
    methods
        function Reset(obj)
        end
        function SavePhaseState(obj)
        end
        function ComputeTotalMobility(obj, ProductionSystem, FluidModel)
            s = ProductionSystem.Reservoir.State.S;
            obj.Mob = FluidModel.ComputePhaseMobilities(s);
            obj.Mobt = sum(obj.Mob,2);
        end
        function UpdateFractionalFlow(obj, ProductionSystem, FluidModel)
            obj.ComputeTotalMobility(ProductionSystem, FluidModel);
            obj.f = obj.Mob(:,1) ./ obj.Mobt;
        end
        function dfdS(obj, ProductionSystem, FluidModel)
            dMob = FluidModel.DMobDS(ProductionSystem.Reservoir.State.S);
            obj.df = (dMob(:,1) .* sum(obj.Mob, 2) - sum(dMob, 2) .* obj.Mob(:,1)) ./ sum(obj.Mob, 2).^2;
        end
        function ComputeTransmissibilities(obj, ProductionSystem, DiscretizationModel)
            % Initialize local variables
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            dx = DiscretizationModel.ReservoirGrid.dx;
            dy = DiscretizationModel.ReservoirGrid.dy;
            dz = DiscretizationModel.ReservoirGrid.dz;
            Ax = DiscretizationModel.ReservoirGrid.Ax;
            Ay = DiscretizationModel.ReservoirGrid.Ay;
            Az = DiscretizationModel.ReservoirGrid.Az;
            K(:, 1) = ProductionSystem.Reservoir.K(:,1) .* obj.Mobt;
            K(:, 2) = ProductionSystem.Reservoir.K(:,2) .* obj.Mobt;
            K(:, 3) = ProductionSystem.Reservoir.K(:,3) .* obj.Mobt;
            
            % Harmonic average of permeability
            kx = reshape(K(:,1), Nx, Ny, Nz);
            ky = reshape(K(:,2), Nx, Ny, Nz);
            kz = reshape(K(:,3), Nx, Ny, Nz);
            Kx = zeros(Nx+1, Ny, Nz);
            Ky = zeros(Nx, Ny+1, Nz);
            Kz = zeros(Nx, Ny, Nz+1);
            Kx(2:Nx,:,:) = 2*kx(1:Nx-1,:,:) .* kx(2:Nx,:,:) ./ (kx(1:Nx-1,:,:) + kx(2:Nx,:,:));
            Ky(:,2:Ny,:) = 2*ky(:,1:Ny-1,:) .* ky(:,2:Ny,:) ./ (ky(:,1:Ny-1,:) + ky(:,2:Ny,:));
            Kz(:,:,2:Nz) = 2*kz(:,:,1:Nz-1) .* kz(:,:,2:Nz) ./ (kz(:,:,1:Nz-1) + kz(:,:,2:Nz));
            
            % Transmissibility
            obj.Tx = zeros(Nx+1, Ny, Nz);
            obj.Ty = zeros(Nx, Ny+1, Nz);
            obj.Tz = zeros(Nx, Ny, Nz+1);
            obj.Tx(2:Nx,:,:) = Ax./dx .* Kx(2:Nx,:,:);
            obj.Ty(:,2:Ny,:) = Ay./dy .* Ky(:,2:Ny,:);
            obj.Tz(:,:,2:Nz) = Az./dz .* Kz(:,:,2:Nz);
        end
        function A = BuildIncompressiblePressureMatrix(obj, DiscretizationModel)
            N =  DiscretizationModel.ReservoirGrid.N;
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            %Construct pressure matrix
            x1 = reshape(obj.Tx(1:Nx,:,:),N,1);
            x2 = reshape(obj.Tx(2:Nx+1,:,:),N,1);
            y1 = reshape(obj.Ty(:,1:Ny,:),N,1);
            y2 = reshape(obj.Ty(:,2:Ny+1,:),N,1);
            z1 = reshape(obj.Tz(:,:,1:Nz),N,1);
            z2 = reshape(obj.Tz(:,:,2:Nz+1),N,1);
            DiagVecs = [-z2, -y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1, -z1];
            DiagIndx = [-Nx*Ny, -Nx,-1,0,1,Nx, Nx*Ny];
            A = spdiags(DiagVecs,DiagIndx,N,N);
        end
        function rhs = BuildIncompressibleRHS(obj, ProductionSystem, DiscretizationModel, FluidModel)
            N = DiscretizationModel.ReservoirGrid.N;
            rhs = zeros(N, 1);
            %Add capillary term to the right-hand side
%             Kw = zeros(2, Grid.Nx, Grid.Ny);
%             Kw(1,:,:) = reshape(Mw, 1, Grid.Nx, Grid.Ny).*K(1,:,:);		% x-direction
%             Kw(2,:,:) = reshape(Mw, 1, Grid.Nx, Grid.Ny).*K(2,:,:);		% y-direction
%             AddPcToPressureSystem(q, S, Fluid, Kw, K(1,:,:), Grid);
        end
        function [A, rhs] = AddWellsToPressureSystem(obj, Wells, K, A, rhs)
            % Injectors
            for i=1:Wells.NofInj                
               [A, rhs] = Wells.Inj(i).AddToPressureSystem(K, A, rhs);
            end
            
            % Producers
            for i=1:Wells.NofProd
                [A, rhs] = Wells.Prod(i).AddToPressureSystem(obj.Mobt, K, A, rhs);
            end
        end
        function ComputeFluxes(obj, ProductionSystem, DiscretizationModel)
            % Initialize local variables 
            p = ProductionSystem.Reservoir.State.p;
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            
            % Compute total fluxes ([m^3/s])
            P = reshape(p, Nx, Ny, Nz);
            obj.U.x = zeros(Nx+1, Ny, Nz);
            obj.U.y = zeros(Nx, Ny+1, Nz);
            obj.U.z = zeros(Nx, Ny, Nz+1);
            obj.U.x(2:Nx,:,:) = (P(1:Nx-1,:,:)-P(2:Nx,:,:)) .* obj.Tx(2:Nx,:,:); %- Ucap.x(2:Nx,:);
            obj.U.y(:,2:Ny,:) = (P(:,1:Ny-1,:)-P(:,2:Ny,:)) .* obj.Ty(:,2:Ny,:); %- Ucap.y(:,2:Ny);
            obj.U.z(:,:,2:Nz) = (P(:,:,1:Nz-1)-P(:,:,2:Nz)) .* obj.Tz(:,:,2:Nz); %- Ucap.y(:,2:Ny);
            
            % Wells total fluxes
            obj.Qwells = ProductionSystem.Wells.TotalFluxes(ProductionSystem.Reservoir, obj.Mobt);
        end
        function conservative = CheckMassConservation(obj, Grid)
            %Checks mass balance in all cells
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            conservative = 1;
            maxUx = max(max(max(obj.U.x)));
            maxUy = max(max(max(obj.U.y)));
            maxUz = max(max(max(obj.U.z)));
            maxU = max([maxUx, maxUy, maxUz]);
            qWells = reshape(obj.Qwells, Nx, Ny, Nz);
            for k=1:Nz
                for j=1:Ny
                    for i=1:Nx
                        Accum = obj.U.x(i,j,k) - obj.U.x(i+1,j,k) + obj.U.y(i,j,k) - obj.U.y(i,j+1,k) + obj.U.z(i,j,k) - obj.U.z(i,j,k+1) + qWells(i,j,k);
                        if (abs(Accum/maxU) > 10^(-5))
                            conservative = 0;
                        end
                    end
                end
            end
        end
        function Residual = BuildTransportResidual(obj, ProductionSystem, DiscretizationModel, dt, State0)
            % Initialise local objects
            pv = ProductionSystem.Reservoir.Por * DiscretizationModel.ReservoirGrid.Volume;
            s = ProductionSystem.Reservoir.State.S;
            s_old = State0.S;
            
            % viscous fluxes matrix
            obj.ViscousMatrix(DiscretizationModel.ReservoirGrid);      
            
            % Compute residual
            Residual = pv/dt * (s - s_old)  - max(obj.Qwells, 0) - obj.V * obj.f;
        end
        function ViscousMatrix(obj, Grid)
            %Builds Upwind Flux matrix
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;                                   
            q = min(obj.Qwells, 0);                        
            % right to left and top to bottom (negative x, y, z)
            Xneg = min(obj.U.x, 0); 
            Yneg = min(obj.U.y, 0);
            Zneg = min(obj.U.z, 0);
            % make them vectors 
            x1 = reshape(Xneg(1:Nx,:,:),N,1);
            y1 = reshape(Yneg(:,1:Ny,:),N,1);
            z1 = reshape(Zneg(:,:,1:Nz),N,1);
            
            % left to right and bottom to top (positive x, y, z)
            Xpos = max(obj.U.x, 0); 
            Ypos = max(obj.U.y, 0); 
            Zpos = max(obj.U.z, 0);
            % make them vectors
            x2 = reshape(Xpos(2:Nx+1,:,:), N, 1);
            y2 = reshape(Ypos(:,2:Ny+1,:), N, 1);
            z2 = reshape(Zpos(:,:,2:Nz+1), N, 1);
            
            % Assemble matrix
            DiagVecs = [z2, y2, x2, q+x1-x2+y1-y2+z1-z2, -x1, -y1, -z1]; % diagonal vectors
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny]; % diagonal index
            obj.V = spdiags(DiagVecs, DiagIndx, N, N);
        end
        function Jacobian = BuildTransportJacobian(obj, ProductionSystem, DiscretizationModel, dt)
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