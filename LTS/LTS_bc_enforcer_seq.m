%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_bc_enforcer_seq < LTS_bc_enforcer
    properties
        ActCells
        BCFluxMatrix
        f
        Vr
        Pc
        Mob2
        ThsBC
        GhsBC
    end
    methods
        function ComputeBoundaryValues(obj, DiscretizationModel, Formulation, CellsSelected)
            
            % store the past values of the fractional flow
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N  = DiscretizationModel.ReservoirGrid.N;
            
            BCFluxes.x = zeros(Nx+1,Ny,Nz);
            BCFluxes.y = zeros(Nx,Ny+1,Nz);
            BCFluxes.z = zeros(Nx,Ny,Nz+1);
            
            ActCellsM = reshape(CellsSelected.ActCells, Nx, Ny, Nz);
            
            % Compute the fluxes values
            Utot = Formulation.Utot;
            
            % right to left and top to bottom (negative x, y, z)
            Xneg = min(Utot.x, 0);
            Yneg = min(Utot.y, 0);
            Zneg = min(Utot.z, 0);
            
            x1 = reshape(Xneg(1:Nx,:,:),N,1);
            y1 = reshape(Yneg(:,1:Ny,:),N,1);
            z1 = reshape(Zneg(:,:,1:Nz),N,1);
            
            % left to right and bottom to top (positive x, y, z)
            Xpos = max(Utot.x, 0);
            Ypos = max(Utot.y, 0);
            Zpos = max(Utot.z, 0);
            
            x2 = reshape(Xpos(2:Nx+1,:,:), N, 1);
            y2 = reshape(Ypos(:,2:Ny+1,:), N, 1);
            z2 = reshape(Zpos(:,:,2:Nz+1), N, 1);
            
            
            %fluxes at the boundaty between an accepted and a refuced cell
            BCFluxes.x(2:Nx,:,:) = (ActCellsM(1:Nx-1,:,:) + ActCellsM(2:Nx,:,:) == 1);
            BCFluxes.y(:,2:Ny,:) = (ActCellsM(:,1:Ny-1,:) + ActCellsM(:,2:Ny,:) == 1);
            BCFluxes.z(:,:,2:Nz) = (ActCellsM(:,:,1:Nz-1) + ActCellsM(:,:,2:Nz) == 1);
            
            
            %select only the bc fluxes
            x1A = reshape(BCFluxes.x(1:Nx,:,:),N,1);
            y1A = reshape(BCFluxes.y(:,1:Ny,:),N,1);
            z1A = reshape(BCFluxes.z(:,:,1:Nz),N,1);
            
            x1 = x1 .* x1A;
            y1 = y1 .* y1A;
            z1 = z1 .* z1A;
            
            %select only the bc fluxes
            x2A = reshape(BCFluxes.x(2:Nx+1,:,:),N,1);
            y2A = reshape(BCFluxes.y(:,2:Ny+1,:),N,1);
            z2A = reshape(BCFluxes.z(:,:,2:Nz+1),N,1);
            
            x2 = x2 .* x2A;
            y2 = y2 .* y2A;
            z2 = z2 .* z2A;
            
            DiagVecs = [z2, y2, x2, +x1-x2+y1-y2+z1-z2, -x1, -y1, -z1]; % diagonal vectors
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            
            % BC for viscosity
            CellsSelected.BCFluxMatrix = spdiags(DiagVecs, DiagIndx, N, N);
            obj.BCFluxMatrix = spdiags(DiagVecs, DiagIndx, N, N);
            
            
            rho =  CellsSelected.rho; 
            UpWind = Formulation.UpWind{1,1};
            
            % Transmissibility matrix construction
            Tx = zeros(Nx+1, Ny, Nz);
            Ty = zeros(Nx, Ny+1, Nz);
            Tz = zeros(Nx, Ny, Nz+1);
            
            % Apply upwind operator
            Mupx = UpWind.x*(obj.Mob2 .* obj.f .* rho);
            Mupy = UpWind.y*(obj.Mob2 .* obj.f .* rho);
            Mupz = UpWind.z*(obj.Mob2 .* obj.f .* rho);
            Mupx = reshape(Mupx, Nx, Ny, Nz);
            Mupy = reshape(Mupy, Nx, Ny, Nz);
            Mupz = reshape(Mupz, Nx, Ny, Nz);
                
            % Transmisibility Matrix
            Tx(2:Nx,:,:)= DiscretizationModel.ReservoirGrid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= DiscretizationModel.ReservoirGrid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= DiscretizationModel.ReservoirGrid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
                        
            % Construct matrix
            x1 = reshape(Tx(1:Nx,:,:), N, 1)   .* x1A;
            x2 = reshape(Tx(2:Nx+1,:,:), N, 1) .* x2A;
            y1 = reshape(Ty(:,1:Ny,:), N, 1)   .* y1A;
            y2 = reshape(Ty(:,2:Ny+1,:), N, 1) .* y2A;
            z1 = reshape(Tz(:,:,1:Nz), N, 1)   .* z1A;
            z2 = reshape(Tz(:,:,2:Nz+1), N, 1) .* z2A;

            
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
            
            %BC for capillary forces
            obj.ThsBC = spdiags(DiagVecs,DiagIndx,N,N);
            
            depth = reshape(DiscretizationModel.ReservoirGrid.Depth, Nx, Ny, Nz);
            
            RhoIntC = Formulation.GravityModel.RhoInt;
            RhoInt.x = (RhoIntC{1,1}.x - RhoIntC{2,1}.x); 
            RhoInt.y = (RhoIntC{1,1}.y - RhoIntC{2,1}.y);
            RhoInt.z = (RhoIntC{1,1}.z - RhoIntC{2,1}.z);
            
            Ugx = zeros(Nx+1, Ny, Nz);
            Ugy = zeros(Nx, Ny+1, Nz);
            Ugz = zeros(Nx, Ny, Nz+1);
            Ugx(2:Nx,:,:) = (depth(1:Nx-1,:,:) - depth(2:Nx,:,:)) .* RhoInt.x(2:Nx,:,:); 
            Ugy(:,2:Ny,:) = (depth(:,1:Ny-1,:) - depth(:,2:Ny,:)) .* RhoInt.y(:,2:Ny,:);      
            Ugz(:,:,2:Nz) = (depth(:,:,1:Nz-1) - depth(:,:,2:Nz)) .* RhoInt.z(:,:,2:Nz);
            
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
            Mupx = A.x*(rho .* obj.Mob2 .* obj.f);
            Mupy = A.y*(rho .* obj.Mob2 .* obj.f);
            Mupz = A.z*(rho .* obj.Mob2 .* obj.f);
            Mupx = reshape(Mupx, Nx, Ny, Nz);
            Mupy = reshape(Mupy, Nx, Ny, Nz);
            Mupz = reshape(Mupz, Nx, Ny, Nz);
            
            % Transmisibility Matrix
            Tx(2:Nx,:,:)= DiscretizationModel.ReservoirGrid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= DiscretizationModel.ReservoirGrid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= DiscretizationModel.ReservoirGrid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
            
            
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
            obj.GhsBC = spdiags(DiagVecs, DiagIndx, N, N);
        end
        
        function ComputeBoundaryValuesSubRef(obj, DiscretizationModel, Formulation, CellsSelected, CellsSelected_old)
            % merge the new accepted fractional flow with those occepted at
            % the previous sub-ref
                        
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N  = DiscretizationModel.ReservoirGrid.N;
            
            CellsSelected.f = Formulation.f .* (CellsSelected_old.ActCells) + ...
                CellsSelected_old.f .* (1 - CellsSelected_old.ActCells);

            CellsSelected.Mob2 = Formulation.Mob(1:N,2) .* (CellsSelected_old.ActCells) + ...
                CellsSelected_old.Mob2 .* (1 - CellsSelected_old.ActCells);

           CellsSelected.Pc = CellsSelected.Pc .* (CellsSelected_old.ActCells) + ...
                CellsSelected_old.Pc .* (1 - CellsSelected_old.ActCells);
            % merge the new fluxes interfaces with those obtained at the
            % previous sub-ref.
            
            BCFluxes.x = zeros(Nx+1,Ny,Nz);
            BCFluxes.y = zeros(Nx,Ny+1,Nz);
            BCFluxes.z = zeros(Nx,Ny,Nz+1);
            ActCellsM = reshape(CellsSelected.ActCells, Nx, Ny, Nz);
            % Compute the fluxes values
            Utot = Formulation.Utot;
            
            % right to left and top to bottom (negative x, y, z)
            Xneg = min(Utot.x, 0);
            Yneg = min(Utot.y, 0);
            Zneg = min(Utot.z, 0);
            
            x1 = reshape(Xneg(1:Nx,:,:),N,1);
            y1 = reshape(Yneg(:,1:Ny,:),N,1);
            z1 = reshape(Zneg(:,:,1:Nz),N,1);
            
            % left to right and bottom to top (positive x, y, z)
            Xpos = max(Utot.x, 0);
            Ypos = max(Utot.y, 0);
            Zpos = max(Utot.z, 0);
            
            x2 = reshape(Xpos(2:Nx+1,:,:), N, 1);
            y2 = reshape(Ypos(:,2:Ny+1,:), N, 1);
            z2 = reshape(Zpos(:,:,2:Nz+1), N, 1);
            
            %fluxes at the boundary between an accepted and a rejected cell
            BCFluxes.x(2:Nx,:,:) = (ActCellsM(1:Nx-1,:,:) + ActCellsM(2:Nx,:,:) == 1);
            BCFluxes.y(:,2:Ny,:) = (ActCellsM(:,1:Ny-1,:) + ActCellsM(:,2:Ny,:) == 1);
            BCFluxes.z(:,:,2:Nz) = (ActCellsM(:,:,1:Nz-1) + ActCellsM(:,:,2:Nz) == 1);
            
            
            %select only the bc fluxes
            x1A = reshape(BCFluxes.x(1:Nx,:,:),N,1);
            y1A = reshape(BCFluxes.y(:,1:Ny,:),N,1);
            z1A = reshape(BCFluxes.z(:,:,1:Nz),N,1);
            
            x1 = x1 .* x1A;
            y1 = y1 .* y1A;
            z1 = z1 .* z1A;
            
            %select only the bc fluxes
            x2A = reshape(BCFluxes.x(2:Nx+1,:,:),N,1);
            y2A = reshape(BCFluxes.y(:,2:Ny+1,:),N,1);
            z2A = reshape(BCFluxes.z(:,:,2:Nz+1),N,1);
            
            x2 = x2 .* x2A;
            y2 = y2 .* y2A;
            z2 = z2 .* z2A;
            
            DiagVecs = [z2, y2, x2, +x1-x2+y1-y2+z1-z2, -x1, -y1, -z1]; % diagonal vectors
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            
            CellsSelected.BCFluxMatrix = spdiags(DiagVecs, DiagIndx, N, N);
            
            rho =  CellsSelected.rho; 
            UpWind = Formulation.UpWind{1,1};
            
            % Transmissibility matrix construction
            Tx = zeros(Nx+1, Ny, Nz);
            Ty = zeros(Nx, Ny+1, Nz);
            Tz = zeros(Nx, Ny, Nz+1);
            
            % Apply upwind operator
            Mupx = UpWind.x*(obj.Mob2 .* obj.f .* rho);
            Mupy = UpWind.y*(obj.Mob2 .* obj.f .* rho);
            Mupz = UpWind.z*(obj.Mob2 .* obj.f .* rho);
            Mupx = reshape(Mupx, Nx, Ny, Nz);
            Mupy = reshape(Mupy, Nx, Ny, Nz);
            Mupz = reshape(Mupz, Nx, Ny, Nz);
                
            % Transmisibility Matrix
            Tx(2:Nx,:,:)= DiscretizationModel.ReservoirGrid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= DiscretizationModel.ReservoirGrid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= DiscretizationModel.ReservoirGrid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
                        
            % Construct matrix
            x1 = reshape(Tx(1:Nx,:,:), N, 1)   .* x1A;
            x2 = reshape(Tx(2:Nx+1,:,:), N, 1) .* x2A;
            y1 = reshape(Ty(:,1:Ny,:), N, 1)   .* y1A;
            y2 = reshape(Ty(:,2:Ny+1,:), N, 1) .* y2A;
            z1 = reshape(Tz(:,:,1:Nz), N, 1)   .* z1A;
            z2 = reshape(Tz(:,:,2:Nz+1), N, 1) .* z2A;

            
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
            
            %BC for capillary forces
            obj.ThsBC = spdiags(DiagVecs,DiagIndx,N,N);
            
            depth = reshape(DiscretizationModel.ReservoirGrid.Depth, Nx, Ny, Nz);
            
            RhoIntC = Formulation.GravityModel.RhoInt;
            RhoInt.x = (RhoIntC{1,1}.x - RhoIntC{2,1}.x); 
            RhoInt.y = (RhoIntC{1,1}.y - RhoIntC{2,1}.y);
            RhoInt.z = (RhoIntC{1,1}.z - RhoIntC{2,1}.z);
            
            Ugx = zeros(Nx+1, Ny, Nz);
            Ugy = zeros(Nx, Ny+1, Nz);
            Ugz = zeros(Nx, Ny, Nz+1);
            Ugx(2:Nx,:,:) = (depth(1:Nx-1,:,:) - depth(2:Nx,:,:)) .* RhoInt.x(2:Nx,:,:); 
            Ugy(:,2:Ny,:) = (depth(:,1:Ny-1,:) - depth(:,2:Ny,:)) .* RhoInt.y(:,2:Ny,:);      
            Ugz(:,:,2:Nz) = (depth(:,:,1:Nz-1) - depth(:,:,2:Nz)) .* RhoInt.z(:,:,2:Nz);
            
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
            Mupx = A.x*(rho .* obj.Mob2 .* obj.f);
            Mupy = A.y*(rho .* obj.Mob2 .* obj.f);
            Mupz = A.z*(rho .* obj.Mob2 .* obj.f);
            Mupx = reshape(Mupx, Nx, Ny, Nz);
            Mupy = reshape(Mupy, Nx, Ny, Nz);
            Mupz = reshape(Mupz, Nx, Ny, Nz);
            
            % Transmisibility Matrix
            Tx(2:Nx,:,:)= DiscretizationModel.ReservoirGrid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= DiscretizationModel.ReservoirGrid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= DiscretizationModel.ReservoirGrid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
            
            
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
            obj.GhsBC = spdiags(DiagVecs, DiagIndx, N, N);
        end
        
        function SetCorrectActiveCells(obj, CellsSelected)
            obj.f = CellsSelected.f;
            obj.Pc = CellsSelected.Pc;
            obj.Mob2 = CellsSelected.Mob2;
            obj.ActCells = CellsSelected.ActCells;
            obj.BCFluxMatrix = CellsSelected.BCFluxMatrix;
        end
        function Residual = AddBC2Residual(obj, Residual, ProductionSystem, Formulation, DiscretizationModel, State0, dt)
            
            N = DiscretizationModel.ReservoirGrid.N; 
            depth = DiscretizationModel.ReservoirGrid.Depth;
            
            rho = ProductionSystem.Reservoir.State.Properties('rho_1').Value;
            % (1) Put to zero the accepted residual cells 
            Residual = (Residual) .* obj.ActCells; 
            % (3) Add the contribution of the accelpted fluxes at the
            % boundary
            Residual = Residual - (obj.BCFluxMatrix * (obj.f .* rho)).* obj.ActCells;
            Residual = Residual - (obj.ThsBC * obj.Pc) .* obj.ActCells;
            Residual = Residual - (obj.GhsBC * depth) .* obj.ActCells;

        end
        function CFL = ComputeCFLNumberLTS(obj, DiscretizationModel, ProductionSystem, dt, Formulation)
            pv = ProductionSystem.Reservoir.Por * DiscretizationModel.ReservoirGrid.Volume;
            
            ActFluxes = Formulation.MatrixAssembler.ActInterfaces;
            
            maxdf = max(abs(Formulation.df .* obj.ActCells));
            maxUx = max(max(max(abs(Formulation.Utot.x .* ActFluxes.x))));
            maxUy = max(max(max(abs(Formulation.Utot.y .* ActFluxes.y))));
            maxUz = max(max(max(abs(Formulation.Utot.z .* ActFluxes.z))));
            CFL = dt * maxdf * (maxUx + maxUy + maxUz) / pv;
        end
        function Jacobian = AddBC2Jacobian(obj, Jacobian, ProductionSystem, Formulation, DiscretizationModel, dt)
            N = DiscretizationModel.ReservoirGrid.N;

            Jacobian = (Jacobian);%.*spdiags(obj.ActCells,0,N,N);
        end
    end
end


 
