%  Cartesian grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 14 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef cartesian_grid < grid_darsim
    properties
        Nx
        Ny
        Nz
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
        pEDFM_alpha_Tx = 0; % Correction for pEDFM Connectivities
        pEDFM_alpha_Ty = 0; % Correction for pEDFM Connectivities
        pEDFM_alpha_Tz = 0; % Correction for pEDFM Connectivities
        condTrans_x % transmisibility of heat conductivity
        condTrans_y % transmisibility of heat conductivity
        condTrans_z % transmisibility of heat conductivity
        I
        J
        K
        Centroids
    end
    methods
        function obj = cartesian_grid(n)
            obj.Nx = n(1);
            obj.Ny = n(2);
            obj.Nz = n(3);
            obj.N = prod(n);
            obj.Active = ones(obj.N, 1);
            obj.ActiveTime = ones(obj.N, 1);
            obj.Tx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            obj.Ty = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            obj.Tz = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            obj.pEDFM_alpha_Tx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            obj.pEDFM_alpha_Ty = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            obj.pEDFM_alpha_Tz = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            obj.condTrans_x = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            obj.condTrans_y = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            obj.condTrans_z = zeros(obj.Nx, obj.Ny, obj.Nz+1);
        end
        function Initialize(obj, Reservoir)
            obj.dx = Reservoir.Length/obj.Nx;
            obj.dy = Reservoir.Width/obj.Ny;
            obj.dz = Reservoir.Thickness/obj.Nz;
            obj.Ax = obj.dy * obj.dz;
            obj.Ay = obj.dx * obj.dz;
            obj.Az = obj.dx * obj.dy;
            obj.Volume = obj.dx * obj.dy * obj.dz;
            obj.AssignNeighbours();
            obj.ComputeRockTransmissibilities(Reservoir.K);
            obj.CoarseFactor = [1, 1, 1];
            obj.DualCoarseType = zeros(obj.N, 1);
            obj.Children = cell(obj.N, 1);
            obj.AddCoordinates();
            obj.Depth = zeros(obj.N, 1);
        end
        function ComputeRockTransmissibilities(obj, K)
            %Harmonic average of permeability.
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
        function AddpEDFMCorrections(obj,pEDFM_alpha_Tx,pEDFM_alpha_Ty,pEDFM_alpha_Tz)
            obj.pEDFM_alpha_Tx = pEDFM_alpha_Tx;
            obj.pEDFM_alpha_Ty = pEDFM_alpha_Ty;
            obj.pEDFM_alpha_Tz = pEDFM_alpha_Tz;
        end
        function CorrectTransmissibilitiesForpEDFM(obj)
            obj.Tx(2:obj.Nx,:,:) = obj.Tx(2:obj.Nx,:,:) .* ( 1 - obj.pEDFM_alpha_Tx(2:obj.Nx,:,:) );
            obj.Ty(:,2:obj.Ny,:) = obj.Ty(:,2:obj.Ny,:) .* ( 1 - obj.pEDFM_alpha_Ty(:,2:obj.Ny,:) );
            obj.Tz(:,:,2:obj.Nz) = obj.Tz(:,:,2:obj.Nz) .* ( 1 - obj.pEDFM_alpha_Tz(:,:,2:obj.Nz) );
%             obj.condTrans_x(2:obj.Nx,:,:) = obj.condTrans_x(2:obj.Nx,:,:) .* ( 1 - obj.pEDFM_alpha_Tx(2:obj.Nx,:,:) );
%             obj.condTrans_y(:,2:obj.Ny,:) = obj.condTrans_y(:,2:obj.Ny,:) .* ( 1 - obj.pEDFM_alpha_Ty(:,2:obj.Ny,:) );
%             obj.condTrans_z(:,:,2:obj.Nz) = obj.condTrans_z(:,:,2:obj.Nz) .* ( 1 - obj.pEDFM_alpha_Tz(:,:,2:obj.Nz) );
        end
        function ComputeHeatConductivitiyTransmissibilities(obj, CondEff)
            % Harmonic average of effective heat conductivity.
            % For now, we assume that heat conductivity is isotropic.
            % This means that CondEff is only one-column vector.
            Cx = reshape(CondEff(:,1), obj.Nx, obj.Ny, obj.Nz);
            Cy = reshape(CondEff(:,1), obj.Nx, obj.Ny, obj.Nz);
            Cz = reshape(CondEff(:,1), obj.Nx, obj.Ny, obj.Nz);
            
            CHx = zeros(obj.Nx+1, obj.Ny, obj.Nz);
            CHy = zeros(obj.Nx, obj.Ny+1, obj.Nz);
            CHz = zeros(obj.Nx, obj.Ny, obj.Nz+1);
            CHx(2:obj.Nx,:,:) = 2*Cx(1:obj.Nx-1,:,:).*Cx(2:obj.Nx,:,:)./(Cx(1:obj.Nx-1,:,:)+Cx(2:obj.Nx,:,:));
            CHy(:,2:obj.Ny,:) = 2*Cy(:,1:obj.Ny-1,:).*Cy(:,2:obj.Ny,:)./(Cy(:,1:obj.Ny-1,:)+Cy(:,2:obj.Ny,:));
            CHz(:,:,2:obj.Nz) = 2*Cz(:,:,1:obj.Nz-1).*Cz(:,:,2:obj.Nz)./(Cy(:,:,1:obj.Nz-1)+Cz(:,:,2:obj.Nz));
            
            %Transmissibility
            obj.condTrans_x(2:obj.Nx,:,:) = obj.Ax./obj.dx.*CHx(2:obj.Nx,:,:);
            obj.condTrans_y(:,2:obj.Ny,:) = obj.Ay./obj.dy.*CHy(:,2:obj.Ny,:);
            obj.condTrans_z(:,:,2:obj.Nz) = obj.Az./obj.dz.*CHz(:,:,2:obj.Nz);
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
            XCM = linspace(obj.dx/2 , obj.Nx*obj.dx - obj.dx/2 , obj.Nx);
            YCM = linspace(obj.dy/2 , obj.Ny*obj.dy - obj.dy/2 , obj.Ny);
            ZCM = linspace(obj.dz/2 , obj.Nz*obj.dz - obj.dz/2 , obj.Nz);
            [x,y,z] = meshgrid(XCM,YCM,ZCM);
            x = permute(x,[2 1 3]);
            y = permute(y,[2 1 3]);
            
            obj.Centroids = [ x(:) , y(:) , z(:) ];
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
        function AssignNeighbours(obj)
            % West Neighbors
            [i,j,k] = meshgrid( 0:obj.Nx-1 , 1:obj.Ny , 1:obj.Nz );  i(i<1) = nan;
            IW = (k-1).*obj.Nx.*obj.Ny + (j-1).*obj.Nx + i;
            IW = reshape(permute(IW,[2 1 3]),obj.N,1);
            
            % East Neighbors
            [i,j,k] = meshgrid( 2:obj.Nx+1 , 1:obj.Ny , 1:obj.Nz );  i(i>obj.Nx) = nan;
            IE = (k-1).*obj.Nx.*obj.Ny + (j-1).*obj.Nx + i;
            IE = reshape(permute(IE,[2 1 3]),obj.N,1);
            
            % South Neighbors
            [i,j,k] = meshgrid( 1:obj.Nx , 0:obj.Ny-1 , 1:obj.Nz );  j(j<1) = nan;
            IS = (k-1).*obj.Nx.*obj.Ny + (j-1).*obj.Nx + i;
            IS = reshape(permute(IS,[2 1 3]),obj.N,1);
            
            % North Neighbors
            [i,j,k] = meshgrid( 1:obj.Nx , 2:obj.Ny+1 , 1:obj.Nz );  j(j>obj.Ny) = nan;
            IN = (k-1).*obj.Nx.*obj.Ny + (j-1).*obj.Nx + i;
            IN = reshape(permute(IN,[2 1 3]),obj.N,1);
            
            % Bottom Neighbors
            [i,j,k] = meshgrid( 1:obj.Nx , 1:obj.Ny , 0:obj.Nz-1 );  k(k<1) = nan;
            IB = (k-1).*obj.Nx.*obj.Ny + (j-1).*obj.Nx + i;
            IB = reshape(permute(IB,[2 1 3]),obj.N,1);
            
            % Top Neighbors
            [i,j,k] = meshgrid( 1:obj.Nx , 1:obj.Ny , 2:obj.Nz+1 );  k(k>obj.Nz) = nan;
            IT = (k-1).*obj.Nx.*obj.Ny + (j-1).*obj.Nx + i;
            IT = reshape(permute(IT,[2 1 3]),obj.N,1);
            
            % Assembling the array for all the neighbors
            obj.Neighbours = [IW,IE,IS,IN,IB,IT];
            obj.Neighbours = num2cell(obj.Neighbours,2);
            obj.Neighbours = cellfun(@(x) x(~isnan(x)),obj.Neighbours ,'UniformOutput' ,false);
        end
    end
end