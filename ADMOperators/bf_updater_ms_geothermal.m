%%  Basis functions updater
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 09 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef bf_updater_ms_geothermal < bf_updater_ms
    properties
    end
    methods
        function ConstructPressureSystem(obj, ProductionSystem, FluidModel, FineGrid, CrossConnections)
            % Builds fine-scale incompressible pressure system.
            % In this function, mobility (Mob) is obtained based on
            % viscosity as no saturation is involved.
            K = ProductionSystem.Reservoir.K;
            mu = ProductionSystem.Reservoir.State.Properties('mu_1').Value;
            Mob = FluidModel.ComputePhaseMobilities(mu);
            obj.A = obj.MediumPressureSystem(FineGrid, K, Mob);
            obj.AddWellsToPressureMatrix(ProductionSystem.Wells, K, Mob, FineGrid.N)
        end
        function ConstructRockTemperatureSystem(obj, ProductionSystem, FineGrid, CrossConnections)
            % Builds fine-scale rock temperature system
            k_cond = ProductionSystem.Reservoir.k_cond;
            obj.A = obj.MediumRockTemperatureSystem(FineGrid(1), k_cond);
        end
        function A_Medium = MediumRockTemperatureSystem(obj, FineGrid, k_cond)
            % Remove high contrast to avoid spikes
            % If MaxContrast is intentially set to 0, avoid this step
            if obj.MaxContrast ~=0
                k_cond_max = max(k_cond(:,1));
                k_cond(k_cond(:,1)./k_cond_max < obj.MaxContrast, 1) = obj.MaxContrast * k_cond_max;
                k_cond(k_cond(:,2)./k_cond_max < obj.MaxContrast, 2) = obj.MaxContrast * k_cond_max;
                k_cond(k_cond(:,3)./k_cond_max < obj.MaxContrast, 3) = obj.MaxContrast * k_cond_max;
            end
            
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
            
            %Harmonic average of conductivity.
            Kx = reshape(k_cond(:,1), Nx, Ny, Nz);
            Ky = reshape(k_cond(:,2), Nx, Ny, Nz);
            Kz = reshape(k_cond(:,3), Nx, Ny, Nz);
            
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
            
            % Correcting for pEDFM connectivities
%             Tx(2:Nx,:,:) = Tx(2:Nx,:,:) .* ( 1 - FineGrid.Tx_Alpha(2:Nx,:,:) );
%             Ty(:,2:Ny,:) = Ty(:,2:Ny,:) .* ( 1 - FineGrid.Ty_Alpha(:,2:Ny,:) );
%             Tz(:,:,2:Nz) = Tz(:,:,2:Nz) .* ( 1 - FineGrid.Tz_Alpha(:,:,2:Nz) );
            
            %Construct rock temperature matrix
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
    end
end