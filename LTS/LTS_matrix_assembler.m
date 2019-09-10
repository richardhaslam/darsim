% Matrix assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_matrix_assembler < matrix_assembler
    properties
        ActInterfaces
    end
    methods
        function ResetActiveInterfaces(obj, DiscretizationModel)
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
                      
            % Reset values to 1
            obj.ActInterfaces.x = ones(Nx+1, Ny, Nz);
            obj.ActInterfaces.y = ones(Nx, Ny+1, Nz);
            obj.ActInterfaces.z = ones(Nx, Ny, Nz+1);
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
            
            % Transmisibility Matrix
            Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
            
            Tx = Tx .* obj.ActInterfaces.x;
            Ty = Ty .* obj.ActInterfaces.y;
            Tz = Tz .* obj.ActInterfaces.z;
            
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
        function [Tpc] = TransmissibilityforPc(obj, Grid, Mob, rho, Pc)
            Pc = reshape(Pc, Nx, Ny, Nz);
            U.x = zeros(Nx+1,Ny,Nz);
            U.y = zeros(Nx,Ny+1,Nz);
            U.z = zeros(Nx,Ny,Nz+1);
            U.x(2:Nx,:,:) = (Pc(1:Nx-1,:,:)-Pc(2:Nx,:,:)) .* Grid.Tx(2:Nx,:,:);
            U.y(:,2:Ny,:) = (Pc(:,1:Ny-1,:)-Pc(:,2:Ny,:)) .* Grid.Ty(:,2:Ny,:);
            U.z(:,:,2:Nz) = (Pc(:,:,1:Nz-1)-Pc(:,:,2:Nz)) .* Grid.Tz(:,:,2:Nz);
            
            %% Use velocity to build upwind operator
            L = reshape((U.x(2:Nx+1,:,:) >= 0), N, 1);
            R = reshape((U.x(1:Nx,:,:) < 0), N, 1);
            B = reshape((U.y(:,2:Ny+1,:) >= 0), N, 1);
            T = reshape((U.y(:,1:Ny,:) < 0), N, 1);
            Down = reshape((U.z(:,:,2:Nz+1) >= 0), N, 1);
            Up = reshape((U.z(:,:,1:Nz) < 0), N, 1);
            
            DiagVecs = [L, R];
            DiagIndx = [0, 1];
            A.x = spdiags(DiagVecs, DiagIndx, N, N);
            DiagVecs = [B, T];
            DiagIndx = [0, Nx];
            A.y = spdiags(DiagVecs, DiagIndx, N, N);
            DiagVecs = [Down, Up];
            DiagIndx = [0, Nx*Ny];
            A.z = spdiags(DiagVecs, DiagIndx, N, N);
            
            % Transmissibility matrix construction
            Tx = zeros(Nx+1, Ny, Nz);
            Ty = zeros(Nx, Ny+1, Nz);
            Tz = zeros(Nx, Ny, Nz+1);
            
            % Apply upwind operator
            Mupx = A.x*(Mob .* rho);
            Mupy = A.y*(Mob .* rho);
            Mupz = A.z*(Mob .* rho);
            Mupx = reshape(Mupx, Nx, Ny, Nz);
            Mupy = reshape(Mupy, Nx, Ny, Nz);
            Mupz = reshape(Mupz, Nx, Ny, Nz);
            
            % Transmisibility Matrix
            Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
            
            Tx = Tx .* obj.ActInterfaces.x;
            Ty = Ty .* obj.ActInterfaces.y;
            Tz = Tz .* obj.ActInterfaces.z;
            
            % Construct matrix
            x1 = reshape(Tx(1:Nx,:,:), N, 1);
            x2 = reshape(Tx(2:Nx+1,:,:), N, 1);
            y1 = reshape(Ty(:,1:Ny,:), N, 1);
            y2 = reshape(Ty(:,2:Ny+1,:), N, 1);
            z1 = reshape(Tz(:,:,1:Nz), N, 1);
            z2 = reshape(Tz(:,:,2:Nz+1), N, 1);
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
            Tpc = spdiags(DiagVecs,DiagIndx,N,N);
        end

        function Ths = TransmissibilityforS(obj, Grid, rho, S)
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            
            S = reshape(S, Nx, Ny, Nz);
            U.x = zeros(Nx+1,Ny,Nz);
            U.y = zeros(Nx,Ny+1,Nz);
            U.z = zeros(Nx,Ny,Nz+1);
            U.x(2:Nx,:,:) = (S(1:Nx-1,:,:) - S(2:Nx,:,:)) .* Grid.Tx(2:Nx,:,:);
            U.y(:,2:Ny,:) = (S(:,1:Ny-1,:) - S(:,2:Ny,:)) .* Grid.Ty(:,2:Ny,:);
            U.z(:,:,2:Nz) = (S(:,:,1:Nz-1) - S(:,:,2:Nz)) .* Grid.Tz(:,:,2:Nz);
            
            %% Use velocity to build upwind operator
            L = reshape((U.x(2:Nx+1,:,:) >= 0), N, 1);
            R = reshape((U.x(1:Nx,:,:) < 0), N, 1);
            B = reshape((U.y(:,2:Ny+1,:) >= 0), N, 1);
            T = reshape((U.y(:,1:Ny,:) < 0), N, 1);
            Down = reshape((U.z(:,:,2:Nz+1) >= 0), N, 1);
            Up = reshape((U.z(:,:,1:Nz) < 0), N, 1);
            
            DiagVecs = [L, R];
            DiagIndx = [0, 1];
            A.x = spdiags(DiagVecs, DiagIndx, N, N);
            DiagVecs = [B, T];
            DiagIndx = [0, Nx];
            A.y = spdiags(DiagVecs, DiagIndx, N, N);
            DiagVecs = [Down, Up];
            DiagIndx = [0, Nx*Ny];
            A.z = spdiags(DiagVecs, DiagIndx, N, N);
            
            % Transmissibility matrix construction
            Tx = zeros(Nx+1, Ny, Nz);
            Ty = zeros(Nx, Ny+1, Nz);
            Tz = zeros(Nx, Ny, Nz+1);
            
            % Apply upwind operator
            Mupx = A.x*(rho);
            Mupy = A.y*(rho);
            Mupz = A.z*(rho);
            Mupx = reshape(Mupx, Nx, Ny, Nz);
            Mupy = reshape(Mupy, Nx, Ny, Nz);
            Mupz = reshape(Mupz, Nx, Ny, Nz);
            
            % Transmisibility Matrix
            Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
            
            Tx = Tx .* obj.ActInterfaces.x;
            Ty = Ty .* obj.ActInterfaces.y;
            Tz = Tz .* obj.ActInterfaces.z;
            
            % Construct matrix
            x1 = reshape(Tx(1:Nx,:,:), N, 1);
            x2 = reshape(Tx(2:Nx+1,:,:), N, 1);
            y1 = reshape(Ty(:,1:Ny,:), N, 1);
            y2 = reshape(Ty(:,2:Ny+1,:), N, 1);
            z1 = reshape(Tz(:,:,1:Nz), N, 1);
            z2 = reshape(Tz(:,:,2:Nz+1), N, 1);
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
            Ths = spdiags(DiagVecs,DiagIndx,N,N);
        end
    end
end

