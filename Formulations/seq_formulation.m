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
        Qwells
        V
        GravityModel
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
            dx = DiscretizationModel.ReservoirGrid.dx;
            dy = DiscretizationModel.ReservoirGrid.dy;
            Ax = DiscretizationModel.ReservoirGrid.Ax;
            Ay = DiscretizationModel.ReservoirGrid.Ay;
            K(:, 1) = ProductionSystem.Reservoir.K(:,1) .* obj.Mobt;
            K(:, 2) = ProductionSystem.Reservoir.K(:,2) .* obj.Mobt;
            
            % Harmonic average of permeability
            kx = reshape(K(:,1), Nx, Ny);
            ky = reshape(K(:,2), Nx, Ny);
            Kx = zeros(Nx+1, Ny);
            Ky = zeros(Nx, Ny+1);
            Kx(2:Nx,:) = 2*kx(1:Nx-1,:) .* kx(2:Nx,:) ./ (kx(1:Nx-1,:) + kx(2:Nx,:));
            Ky(:,2:Ny) = 2*ky(:,1:Ny-1) .* ky(:,2:Ny) ./ (ky(:,1:Ny-1) + ky(:,2:Ny));
            
            % Transmissibility
            obj.Tx = zeros(Nx+1, Ny);
            obj.Ty = zeros(Nx, Ny+1);
            obj.Tx(2:Nx,:) = Ax./dx .* Kx(2:Nx,:);
            obj.Ty(:,2:Ny) = Ay./dy .* Ky(:,2:Ny);
        end
        function A = BuildIncompressiblePressureMatrix(obj, DiscretizationModel)
            N =  DiscretizationModel.ReservoirGrid.N;
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            %Construct pressure matrix
            x1 = reshape(obj.Tx(1:Nx,:),N,1);
            x2 = reshape(obj.Tx(2:Nx+1,:),N,1);
            y1 = reshape(obj.Ty(:,1:Ny),N,1);
            y2 = reshape(obj.Ty(:,2:Ny+1),N,1);
            DiagVecs = [-y2,-x2,y2+x2+y1+x1,-x1,-y1];
            DiagIndx = [-Nx,-1,0,1,Nx];
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
            
            % Compute total fluxes ([m^3/s])
            P = reshape(p, Nx, Ny, 1);
            obj.U.x = zeros(Nx+1,Ny, 1);
            obj.U.y = zeros(Nx,Ny+1, 1);
            obj.U.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)) .* obj.Tx(2:Nx,:); %- Ucap.x(2:Nx,:);
            obj.U.y(:,2:Ny) = (P(:,1:Ny-1)-P(:,2:Ny)) .* obj.Ty(:,2:Ny); %- Ucap.y(:,2:Ny);
            
            % Wells total fluxes
            obj.Qwells = ProductionSystem.Wells.TotalFluxes(ProductionSystem.Reservoir, obj.Mobt);
        end
        function conservative = CheckMassConservation(obj, Grid)
            %Checks mass balance in all cells
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            conservative = 1;
            maxUx = max(max(obj.U.x));
            maxUy = max(max(obj.U.y));
            maxU = max(maxUx, maxUy);
            qWells = reshape(obj.Qwells, Nx, Ny);
            for i=1:Nx
                for j=1:Ny
                    Accum = obj.U.x(i,j) - obj.U.x(i+1,j) + obj.U.y(i,j) - obj.U.y(i,j+1) + qWells(i,j);
                    if (abs(Accum/maxU) > 10^(-5))
                        conservative = 0;
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
            Residual = pv/dt * (s - s_old)  - max(obj.Qwells,0) - obj.V*obj.f;
        end
        function ViscousMatrix(obj, Grid)
            %Builds Upwind Flux matrix
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            N = Grid.N;                                    % number of cells
            q = min(obj.Qwells,0);                        % production >> it's a negative flux
            % right to left and top to bottom (negative x and y)
            XN = min(obj.U.x,0); 
            YN = min(obj.U.y,0); 
            % make them vectors 
            x1 = reshape(XN(1:Nx,:),N,1);
            y1 = reshape(YN(:,1:Ny),N,1);
            
            % left to right and bottom to top (positive x and y)
            XP = max(obj.U.x,0); 
            YP = max(obj.U.y,0); 
            % make them vectors
            x2 = reshape(XP(2:Nx+1,:),N,1);
            y2 = reshape(YP(:,2:Ny+1),N,1);
            
            % Assemble matrix
            DiagVecs = [y2, x2, q+x1-x2+y1-y2, -x1, -y1]; % diagonal vectors
            DiagIndx = [-Nx, -1, 0, 1, Nx]; % diagonal index
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