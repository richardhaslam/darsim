% Matrix assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef matrix_assembler < handle
    properties
        
    end
    methods
        function [Tph, Gph] = TransmissibilityMatrix(obj, Grid, Upwind, Mob, rho, RhoInt)
            switch class(Grid)
                case('corner_point_grid')
                    N_Faces = length(Grid.Trans);
                    
                    %% Mass Flux Transmissibility and Gravity Matrix
                    % Apply Upwind Operator
                    Mob_rho_Upwind = Upwind * (Mob.*rho);
                    
                    % Transmisibility Matrix
                    Tph = Grid.ConnectivityMatrix * spdiags(Grid.Trans.*Mob_rho_Upwind,0,N_Faces,N_Faces) * Grid.ConnectivityMatrix';
                    
                    % Gravity Matrix
                    Gph = Grid.ConnectivityMatrix * spdiags(Grid.Trans.*RhoInt        ,0,N_Faces,N_Faces) * Grid.ConnectivityMatrix';

                case('cartesian_grid')
                    Nx = Grid.Nx;
                    Ny = Grid.Ny;
                    Nz = Grid.Nz;
                    N = Grid.N;
                    
                    %% Mass Flux Transmissibility and Gravity Matrix
                    Tx = zeros(Nx+1, Ny, Nz);
                    Ty = zeros(Nx, Ny+1, Nz);
                    Tz = zeros(Nx, Ny, Nz+1);
                    
                    % Apply Upwind Operator
                    Mupx = Upwind.x*(Mob .* rho);
                    Mupy = Upwind.y*(Mob .* rho);
                    Mupz = Upwind.z*(Mob .* rho);
                    Mupx = reshape(Mupx, Nx, Ny, Nz);
                    Mupy = reshape(Mupy, Nx, Ny, Nz);
                    Mupz = reshape(Mupz, Nx, Ny, Nz);
                    
                    % Transmisibility Matrix
                    Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:) .* Mupx(1:Nx-1,:,:);
                    Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:) .* Mupy(:,1:Ny-1,:);
                    Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz) .* Mupz(:,:,1:Nz-1);
                    Tph = obj.ReshapeCartesianTransmisibility(Nx, Ny, Nz, N, Tx, Ty, Tz);
                    
                    % Gravity Matrix
                    Tx(2:Grid.Nx,:,:)= Tx(2:Grid.Nx,:,:) .* RhoInt.x(2:Grid.Nx,:,:);
                    Ty(:,2:Grid.Ny,:)= Ty(:,2:Grid.Ny,:) .* RhoInt.y(:,2:Grid.Ny,:);
                    Tz(:,:,2:Grid.Nz)= Tz(:,:,2:Grid.Nz) .* RhoInt.z(:,:,2:Grid.Nz);
                    Gph = obj.ReshapeCartesianTransmisibility(Nx, Ny, Nz, N, Tx, Ty, Tz);
            end
        end
        function T_Matrix = ReshapeCartesianTransmisibility(obj, Nx, Ny, Nz, N, Tx, Ty, Tz) % remove function - end
                x1 = reshape( Tx(1:Nx  ,:     ,:     ) , N , 1 );
                x2 = reshape( Tx(2:Nx+1,:     ,:     ) , N , 1 );
                y1 = reshape( Ty(:     ,1:Ny  ,:     ) , N , 1 );
                y2 = reshape( Ty(:     ,2:Ny+1,:     ) , N , 1 );
                z1 = reshape( Tz(:     ,:     ,1:Nz  ) , N , 1 );
                z2 = reshape( Tz(:     ,:     ,2:Nz+1) , N , 1 );
                DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
                DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
                T_Matrix = spdiags(DiagVecs,DiagIndx,N,N);
            end
        function Vr = ViscousMatrix(obj, Grid, Qwells, Utot)
            %Builds Upwind Flux matrix
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            q = min(Qwells, 0);
            
            % right to left and top to bottom (negative x, y, z)
            Xneg = min(Utot.x, 0);
            Yneg = min(Utot.y, 0);
            Zneg = min(Utot.z, 0);
            
            % make them vectors
            x1 = reshape(Xneg(1:Nx,:,:),N,1);
            y1 = reshape(Yneg(:,1:Ny,:),N,1);
            z1 = reshape(Zneg(:,:,1:Nz),N,1);
            
            % left to right and bottom to top (positive x, y, z)
            Xpos = max(Utot.x, 0);
            Ypos = max(Utot.y, 0);
            Zpos = max(Utot.z, 0);
            % make them vectors
            x2 = reshape(Xpos(2:Nx+1,:,:), N, 1);
            y2 = reshape(Ypos(:,2:Ny+1,:), N, 1);
            z2 = reshape(Zpos(:,:,2:Nz+1), N, 1);
            
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
            

            % Transmissibility matrix construction
            Tx = zeros(Nx+1, Ny, Nz);
            Ty = zeros(Nx, Ny+1, Nz);
            Tz = zeros(Nx, Ny, Nz+1);
            
    
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
 
             % Gravity Matrix
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
            x1 = reshape(Tx(1:Nx,:,:), N, 1);
            x2 = reshape(Tx(2:Nx+1,:,:), N, 1);
            y1 = reshape(Ty(:,1:Ny,:), N, 1);
            y2 = reshape(Ty(:,2:Ny+1,:), N, 1);           
            z1 = reshape(Tz(:,:,1:Nz), N, 1);
            z2 = reshape(Tz(:,:,2:Nz+1), N, 1);
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
            Ghs = spdiags(DiagVecs, DiagIndx, N, N);
        end
    end
end

