%  Cartesian grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author:
%TU Delft
%Created:
%Last modified:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef corner_point_grid < grid_darsim
    properties
        Nx
        Ny
        Nz
        N_ActiveCells
        N_InternalFaces
        N_ExternalFaces
        dx
        dy
        dz
        Ax
        Ay
        Az
        Volume
        Tx
        Ty
        Tz
        Tx_Alpha = 0; % Correction for pEDFM Connectivities
        Ty_Alpha = 0; % Correction for pEDFM Connectivities
        Tz_Alpha = 0; % Correction for pEDFM Connectivities
        THx % transmisibility of heat conductivity
        THy % transmisibility of heat conductivity
        THz % transmisibility of heat conductivity
        I
        J
        K
    end
    methods
        function obj = corner_point_grid(n)
            obj.Nx = n(1);
            obj.Ny = n(2);
            obj.Nz = n(3);
            %obj.N_ActiveCells = ...;
            obj.N = prod(n);
            obj.Active = ones(obj.N, 1);
            obj.ActiveTime = ones(obj.N, 1);
            obj.Tx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            obj.Ty = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            obj.Tz = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            obj.THx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            obj.THy = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            obj.THz = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            obj.Tx_Alpha = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            obj.Ty_Alpha = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            obj.Tz_Alpha = zeros(obj.Nx, obj.Ny, obj.Nz+1);
        end
        function Initialize(obj, Reservoir)
            obj.dx = Reservoir.Length/obj.Nx;
            obj.dy = Reservoir.Width/obj.Ny;
            obj.dz = Reservoir.Thickness/obj.Nz;
            obj.Ax = obj.dy * obj.dz;
            obj.Ay = obj.dx * obj.dz;
            obj.Az = obj.dx * obj.dy;
            obj.Volume = obj.dx * obj.dy * obj.dz;
            obj.ComputeRockTransmissibilities(Reservoir.K);
%             obj.ComputeRockHeatConductivities(Reservoir.K_Cond_eff);
            obj.CoarseFactor = [1, 1, 1];
            obj.Children = cell(obj.N, 1);
            obj.GrandChildren = cell(obj.N, 1);
            obj.AddCoordinates();
            obj.Depth = zeros(obj.N, 1);
        end
        function AddpEDFMCorrections(obj,Tx_Alpha,Ty_Alpha,Tz_Alpha)
            obj.Tx_Alpha = Tx_Alpha;
            obj.Ty_Alpha = Ty_Alpha;
            obj.Tz_Alpha = Tz_Alpha;
        end
        function ComputeRockTransmissibilities(obj, K)
            %Harmonic average of permeability.
            obj.Trans = zeros(obj.N_InternalFaces);
            for i = 1 : obj.N_InternalFaces
                %obj.Trans(i) = ...;
            end
            Kx = reshape(K(:,1), obj.Nx, obj.Ny, obj.Nz);
            Ky = reshape(K(:,2), obj.Nx, obj.Ny, obj.Nz);
            Kz = reshape(K(:,3), obj.Nx, obj.Ny, obj.Nz);
            
            KHx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            KHy = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            KHz = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            KHx(2:obj.Nx,:,:) = 2*Kx(1:obj.Nx-1,:,:).*Kx(2:obj.Nx,:,:)./(Kx(1:obj.Nx-1,:,:)+Kx(2:obj.Nx,:,:));
            KHy(:,2:obj.Ny,:) = 2*Ky(:,1:obj.Ny-1,:).*Ky(:,2:obj.Ny,:)./(Ky(:,1:obj.Ny-1,:)+Ky(:,2:obj.Ny,:));
            KHz(:,:,2:obj.Nz) = 2*Kz(:,:,1:obj.Nz-1).*Kz(:,:,2:obj.Nz)./(Ky(:,:,1:obj.Nz-1)+Kz(:,:,2:obj.Nz));
            
            %Transmissibility
            obj.Tx(2:obj.Nx,:,:) = obj.Ax./obj.dx.*KHx(2:obj.Nx,:,:);
            obj.Ty(:,2:obj.Ny,:) = obj.Ay./obj.dy.*KHy(:,2:obj.Ny,:);
            obj.Tz(:,:,2:obj.Nz) = obj.Az./obj.dz.*KHz(:,:,2:obj.Nz);
        end
        function CorrectTransmissibilitiesForpEDFM(obj)
            obj.Tx(2:obj.Nx,:,:) = obj.Tx(2:obj.Nx,:,:) .* ( 1 - obj.Tx_Alpha(2:obj.Nx,:,:) );
            obj.Ty(:,2:obj.Ny,:) = obj.Ty(:,2:obj.Ny,:) .* ( 1 - obj.Ty_Alpha(:,2:obj.Ny,:) );
            obj.Tz(:,:,2:obj.Nz) = obj.Tz(:,:,2:obj.Nz) .* ( 1 - obj.Tz_Alpha(:,:,2:obj.Nz) );
            
            obj.THx(2:obj.Nx,:,:) = obj.THx(2:obj.Nx,:,:) .* ( 1 - obj.Tx_Alpha(2:obj.Nx,:,:) );
            obj.THy(:,2:obj.Ny,:) = obj.THy(:,2:obj.Ny,:) .* ( 1 - obj.Ty_Alpha(:,2:obj.Ny,:) );
            obj.THz(:,:,2:obj.Nz) = obj.THz(:,:,2:obj.Nz) .* ( 1 - obj.Tz_Alpha(:,:,2:obj.Nz) );
        end
        function ComputeRockHeatConductivities(obj, K_Cond)
            %Harmonic average of permeability.
            Kx = reshape(K_Cond(:,1), obj.Nx, obj.Ny, obj.Nz);
            Ky = reshape(K_Cond(:,2), obj.Nx, obj.Ny, obj.Nz);
            Kz = reshape(K_Cond(:,3), obj.Nx, obj.Ny, obj.Nz);
            
            KHx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            KHy = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            KHz = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            KHx(2:obj.Nx,:,:) = 2*Kx(1:obj.Nx-1,:,:).*Kx(2:obj.Nx,:,:)./(Kx(1:obj.Nx-1,:,:)+Kx(2:obj.Nx,:,:));
            KHy(:,2:obj.Ny,:) = 2*Ky(:,1:obj.Ny-1,:).*Ky(:,2:obj.Ny,:)./(Ky(:,1:obj.Ny-1,:)+Ky(:,2:obj.Ny,:));
            KHz(:,:,2:obj.Nz) = 2*Kz(:,:,1:obj.Nz-1).*Kz(:,:,2:obj.Nz)./(Ky(:,:,1:obj.Nz-1)+Kz(:,:,2:obj.Nz));
            
            %Transmissibility
            obj.THx(2:obj.Nx,:,:) = obj.Ax./obj.dx.*KHx(2:obj.Nx,:,:);
            obj.THy(:,2:obj.Ny,:) = obj.Ay./obj.dy.*KHy(:,2:obj.Ny,:);
            obj.THz(:,:,2:obj.Nz) = obj.Az./obj.dz.*KHz(:,:,2:obj.Nz);
        end
        function AddCoordinates(obj)
            % Add I, J coordinates to Grid
            obj.I = ones(obj.N, 1);
            obj.J = ones(obj.N, 1);
            obj.K = ones(obj.N, 1);
            Jindexes = 1:1:obj.Ny;
            for k=1:obj.Nz
                for i=1:obj.Ny
                    a = obj.Nx*(i-1)+ (k-1)*obj.Nx*obj.Ny + 1;
                    obj.I(a:a+obj.Nx-1) = 1:1:obj.Nx;
                    obj.J(a:a+obj.Nx-1) = Jindexes(i)*ones(obj.Nx, 1);
                end
            end
            for i=1:obj.Nz
                obj.K((i-1)*obj.Nx*obj.Ny+1:i*obj.Nx*obj.Ny) = i*ones(obj.Nx*obj.Ny, 1);
            end
        end
        function AddGridCoordinates(obj)
            Xi = linspace(0,obj.dx*obj.Nx,obj.Nx+1)';
            Yi = linspace(0,obj.dy*obj.Ny,obj.Ny+1)';
            Zi = linspace(0,obj.dz*obj.Nz,obj.Nz+1)';
            [x,y,z] = meshgrid(Xi,Yi,Zi);
            x = permute(x,[2 1 3]);
            y = permute(y,[2 1 3]);
            z = permute(z,[2 1 3]);
            obj.GridCoords = [x(:),y(:),z(:)];
        end
        function ComputeDepth(obj, alpha, Thickness)
            x_centres = (obj.I - 1/2) * obj.dx;
            y_centres = (obj.J - 1/2) * obj.dy;
            z_centres = (obj.K - 1/2) * obj.dz;
            obj.Depth = Thickness - z_centres;
        end
    end
end