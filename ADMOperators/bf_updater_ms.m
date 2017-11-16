%  Basis functions updater
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 09 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef bf_updater_ms < bf_updater
    properties
        MaxContrast
    end
    methods
        function ConstructPressureSystem(obj, ProductionSystem, FluidModel, FineGrid, CrossConnections)
            % Builds fine-scale incompressible pressure system
            K = ProductionSystem.Reservoir.K;
            Sm = ProductionSystem.Reservoir.State.Properties('S_1').Value;
            Mob = FluidModel.ComputePhaseMobilities(Sm);
            obj.A = obj.MediumPressureSystem(FineGrid, K, Mob);
            obj.AddWellsToPressureMatrix(ProductionSystem.Wells, K, Mob, FineGrid.N)
        end
        function A_Medium = MediumPressureSystem(obj, FineGrid, K, Mob)
            % Remove high contrast to avoid spikes
            lambdaMax = max(K(:,1));
            K(K(:,1)./lambdaMax < obj.MaxContrast, 1) = obj.MaxContrast * lambdaMax;
            K(K(:,2)./lambdaMax < obj.MaxContrast, 2) = obj.MaxContrast * lambdaMax;
            K(K(:,3)./lambdaMax < obj.MaxContrast, 3) = obj.MaxContrast * lambdaMax;
            
            % Initialize local variables
            Nx = FineGrid.Nx;
            Ny = FineGrid.Ny;
            Nz = FineGrid.Nz;
            N = FineGrid.N;
            dx = FineGrid.dx;
            dy = FineGrid.dy;
            dz = FineGrid.dz;
            Ax = FineGrid.Ax;
            Ay = FineGrid.Ay;
            Az = FineGrid.Az;
            Mobt = sum(Mob,2);
            
            K(:, 1) = K(:,1) .* Mobt;
            K(:, 2) = K(:,2) .* Mobt;
            K(:, 3) = K(:,3) .* Mobt;
            
            %Harmonic average of permeability.
            Kx = reshape(K(:,1), Nx, Ny, Nz);
            Ky = reshape(K(:,2), Nx, Ny, Nz);
            Kz = reshape(K(:,3), Nx, Ny, Nz);
            
            KHx = zeros(Nx+1, Ny, Nz);
            KHy = zeros(Nx, Ny+1, Nz);
            KHz = zeros(Nx, Ny, Nz+1);
            KHx(2:Nx,:,:) = 2*Kx(1:Nx-1,:,:).*Kx(2:Nx,:,:)./(Kx(1:Nx-1,:,:)+Kx(2:Nx,:,:));
            KHy(:,2:Ny,:) = 2*Ky(:,1:Ny-1,:).*Ky(:,2:Ny,:)./(Ky(:,1:Ny-1,:)+Ky(:,2:Ny,:));
            KHz(:,:,2:Nz) = 2*Kz(:,:,1:Nz-1).*Kz(:,:,2:Nz)./(Ky(:,:,1:Nz-1)+Kz(:,:,2:Nz));
            
            %Transmissibility
            Tx = zeros(Nx+1, Ny, Nz);
            Ty = zeros(Nx, Ny+1, Nz);
            Tz = zeros(Nx, Ny, Nz+1);
            Tx(2:Nx,:,:) = Ax./dx.*KHx(2:Nx,:,:);
            Ty(:,2:Ny,:) = Ay./dy.*KHy(:,2:Ny,:);
            Tz(:,:,2:Nz) = Az./dz.*KHz(:,:,2:Nz);
            
            %Construct pressure matrix
            x1 = reshape(Tx(1:Nx,:,:),N,1);
            x2 = reshape(Tx(2:Nx+1,:,:),N,1);
            y1 = reshape(Ty(:,1:Ny,:),N,1);
            y2 = reshape(Ty(:,2:Ny+1,:),N,1);
            z1 = reshape(Tz(:,:,1:Nz),N,1);
            z2 = reshape(Tz(:,:,2:Nz+1),N,1); 
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
            A_Medium = spdiags(DiagVecs, DiagIndx, N, N);
        end
        function AddWellsToPressureMatrix(obj, Wells, K, Mob, N)
            %% Add Wells in residual form
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            dq = zeros(N, 1);
            
            %Injectors
            for i=1:length(Inj)
                c = Inj(i).Cells;
                dq(c) = Inj(i).PI * sum(Inj(i).Mob, 2) * K(c, 1);
            end
            
            %Producers
            Mobt = sum(Mob, 2);
            for i=1:length(Prod)
                c = Prod(i).Cells;
                dq(c) = Prod(i).PI * Mobt(c, :) .* K(c,1);
            end
 
            W = spdiags(dq, 0, N, N);
            obj.A = obj.A + W;
        end
    end
end