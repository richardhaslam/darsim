% Matrix assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_matrix_assembler < matrix_assembler
    properties
        ActInterfaces
        ActCells
    end
    methods
        function ResetActiveInterfaces(obj, DiscretizationModel)
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N = DiscretizationModel.ReservoirGrid.N;
                      
            % Reset values to 1
            obj.ActInterfaces.x = ones(Nx+1, Ny, Nz);
            obj.ActInterfaces.y = ones(Nx, Ny+1, Nz);
            obj.ActInterfaces.z = ones(Nx, Ny, Nz+1);

            obj.ActCells = ones(N,1);
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
        function Vr = ViscousMatrix(obj, Grid, Qwells, Utot)
            %Builds Upwind Flux matrix
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;                                   
            q = min(Qwells, 0) .* obj.ActCells;  
             
            % right to left and top to bottom (negative x, y, z)
            Xneg = min(Utot.x, 0); 
            Yneg = min(Utot.y, 0);
            Zneg = min(Utot.z, 0);
            
            % make them vectors 
            x1 = reshape(Xneg(1:Nx,:,:),N,1);
            y1 = reshape(Yneg(:,1:Ny,:),N,1);
            z1 = reshape(Zneg(:,:,1:Nz),N,1);
            
            %select only the active fluxes
            x1A = reshape(obj.ActInterfaces.x(1:Nx,:,:),N,1);
            y1A = reshape(obj.ActInterfaces.y(:,1:Ny,:),N,1);
            z1A = reshape(obj.ActInterfaces.z(:,:,1:Nz),N,1);
            
            x1 = x1 .* x1A;
            y1 = y1 .* y1A;
            z1 = z1 .* z1A;
            
            % left to right and bottom to top (positive x, y, z)
            Xpos = max(Utot.x, 0); 
            Ypos = max(Utot.y, 0); 
            Zpos = max(Utot.z, 0);
            % make them vectors
            x2 = reshape(Xpos(2:Nx+1,:,:), N, 1);
            y2 = reshape(Ypos(:,2:Ny+1,:), N, 1);
            z2 = reshape(Zpos(:,:,2:Nz+1), N, 1);
            
            %select only the active fluxes
            x2A = reshape(obj.ActInterfaces.x(2:Nx+1,:,:),N,1);
            y2A = reshape(obj.ActInterfaces.y(:,2:Ny+1,:),N,1);
            z2A = reshape(obj.ActInterfaces.z(:,:,2:Nz+1),N,1);
            
            x2 = x2 .* x2A;
            y2 = y2 .* y2A;
            z2 = z2 .* z2A;
            
            % Assemble matrix
            DiagVecs = [z2, y2, x2, +x1-x2+y1-y2+z1-z2+q, -x1, -y1, -z1]; % diagonal vectors
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny]; % diagonal index
            Vr = spdiags(DiagVecs, DiagIndx, N, N);
        end
        function [Ths, UPc, Ghs] = TransmissibilityforS(obj, Grid, UpWind, vect, rho, Pc, RhoIntC)
            
            RhoInt.x = (RhoIntC{1,1}.x - RhoIntC{2,1}.x); 
            RhoInt.y = (RhoIntC{1,1}.y - RhoIntC{2,1}.y);
            RhoInt.z = (RhoIntC{1,1}.z - RhoIntC{2,1}.z);
            
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            
            Pc = reshape(Pc, Nx, Ny, Nz);
            depth = reshape(Grid.Depth, Nx, Ny, Nz);
            
            %
            Ugx = zeros(Nx+1, Ny, Nz);
            Ugy = zeros(Nx, Ny+1, Nz);
            Ugz = zeros(Nx, Ny, Nz+1);
            Ugx(2:Nx,:,:) = (depth(1:Nx-1,:,:) - depth(2:Nx,:,:)) .* RhoInt.x(2:Nx,:,:); 
            Ugy(:,2:Ny,:) = (depth(:,1:Ny-1,:) - depth(:,2:Ny,:)) .* RhoInt.y(:,2:Ny,:);      
            Ugz(:,:,2:Nz) = (depth(:,:,1:Nz-1) - depth(:,:,2:Nz)) .* RhoInt.z(:,:,2:Nz);
            
            %% Compute  fluxes ([m^3/s])
            UPc.x = zeros(Nx+1,Ny,Nz);
            UPc.y = zeros(Nx,Ny+1,Nz);
            UPc.z = zeros(Nx,Ny,Nz+1);
            UPc.x(2:Nx,:,:) = -(Pc(1:Nx-1,:,:)-Pc(2:Nx,:,:)) .* Grid.Tx(2:Nx,:,:) - Grid.Tx(2:Nx,:,:) .* Ugx(2:Nx,:,:);
            UPc.y(:,2:Ny,:) = -(Pc(:,1:Ny-1,:)-Pc(:,2:Ny,:)) .* Grid.Ty(:,2:Ny,:) - Grid.Ty(:,2:Ny,:) .* Ugy(:,2:Ny,:);
            UPc.z(:,:,2:Nz) = -(Pc(:,:,1:Nz-1)-Pc(:,:,2:Nz)) .* Grid.Tz(:,:,2:Nz) - Grid.Tz(:,:,2:Nz) .* Ugz(:,:,2:Nz);

            
            Pc = reshape(Pc, Nx, Ny, Nz);
            %% Compute  fluxes ([m^3/s])
            UPc.x = zeros(Nx+1,Ny,Nz);
            UPc.y = zeros(Nx,Ny+1,Nz);
            UPc.z = zeros(Nx,Ny,Nz+1);
            UPc.x(2:Nx,:,:) = -(Pc(1:Nx-1,:,:)-Pc(2:Nx,:,:)) .* Grid.Tx(2:Nx,:,:);
            UPc.y(:,2:Ny,:) = -(Pc(:,1:Ny-1,:)-Pc(:,2:Ny,:)) .* Grid.Ty(:,2:Ny,:);
            UPc.z(:,:,2:Nz) = -(Pc(:,:,1:Nz-1)-Pc(:,:,2:Nz)) .* Grid.Tz(:,:,2:Nz);
            
 
            UPc.x = UPc.x .* obj.ActInterfaces.x;
            UPc.y = UPc.y .* obj.ActInterfaces.y;
            UPc.z = UPc.z .* obj.ActInterfaces.z;
            % Transmissibility matrix construction
            Tx = zeros(Nx+1, Ny, Nz);
            Ty = zeros(Nx, Ny+1, Nz);
            Tz = zeros(Nx, Ny, Nz+1);
            
            % Apply upwind operator
            Mupx = UpWind.x*(vect .* rho);
            Mupy = UpWind.y*(vect .* rho);
            Mupz = UpWind.z*(vect .* rho);
            Mupx = reshape(Mupx, Nx, Ny, Nz);
            Mupy = reshape(Mupy, Nx, Ny, Nz);
            Mupz = reshape(Mupz, Nx, Ny, Nz);
            
            % Transmisibility Matrix
            Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);


            x1A = reshape(obj.ActInterfaces.x(1:Nx,:,:),N,1);
            y1A = reshape(obj.ActInterfaces.y(:,1:Ny,:),N,1);
            z1A = reshape(obj.ActInterfaces.z(:,:,1:Nz),N,1);
            
            x2A = reshape(obj.ActInterfaces.x(2:Nx+1,:,:),N,1);
            y2A = reshape(obj.ActInterfaces.y(:,2:Ny+1,:),N,1);
            z2A = reshape(obj.ActInterfaces.z(:,:,2:Nz+1),N,1);
            
            % Construct matrix
            x1 = reshape(Tx(1:Nx,:,:), N, 1)  .* x1A;
            x2 = reshape(Tx(2:Nx+1,:,:), N, 1).* x2A;
            y1 = reshape(Ty(:,1:Ny,:), N, 1)  .* y1A;
            y2 = reshape(Ty(:,2:Ny+1,:), N, 1).* y2A;
            z1 = reshape(Tz(:,:,1:Nz), N, 1)  .* z1A;
            z2 = reshape(Tz(:,:,2:Nz+1), N, 1).* z2A;
            
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
            Ths = spdiags(DiagVecs,DiagIndx,N,N);
            
            
            %% Use velocity to build upwind operator
       
            L = reshape((Ugx(2:Nx+1,:,:) >= 0), N, 1);
            R = reshape((Ugx(1:Nx,:,:) < 0), N, 1);
            B = reshape((Ugy(:,2:Ny+1,:) >= 0), N, 1);
            T = reshape((Ugy(:,1:Ny,:) < 0), N, 1);
            Down = reshape((Ugz(:,:,2:Nz+1) >= 0), N, 1);
            Up = reshape((Ugz(:,:,1:Nz) < 0), N, 1);
            
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
            Mupx = A.x*(rho .* vect);
            Mupy = A.y*(rho .* vect);
            Mupz = A.z*(rho .* vect);
            Mupx = reshape(Mupx, Nx, Ny, Nz);
            Mupy = reshape(Mupy, Nx, Ny, Nz);
            Mupz = reshape(Mupz, Nx, Ny, Nz);
            
            % Transmisibility Matrix
            Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
            
            
            Tx(2:Nx,:,:)= Tx(2:Nx,:,:) .* RhoInt.x(2:Nx,:,:);
            Ty(:,2:Ny,:)= Ty(:,2:Ny,:) .* RhoInt.y(:,2:Ny,:);
            Tz(:,:,2:Nz)= Tz(:,:,2:Nz) .* RhoInt.z(:,:,2:Nz);
            % Construct matrix
            x1 = reshape(Tx(1:Nx,:,:), N, 1)   .* x1A;
            x2 = reshape(Tx(2:Nx+1,:,:), N, 1) .* x2A;
            y1 = reshape(Ty(:,1:Ny,:), N, 1)   .* y1A;
            y2 = reshape(Ty(:,2:Ny+1,:), N, 1) .* y2A;           
            z1 = reshape(Tz(:,:,1:Nz), N, 1)   .* z1A;
            z2 = reshape(Tz(:,:,2:Nz+1), N, 1) .* z2A;
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
            Ghs = spdiags(DiagVecs, DiagIndx, N, N);
  
        end
    end
end

