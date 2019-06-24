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
            
            CellsSelected.BCFluxMatrix = spdiags(DiagVecs, DiagIndx, N, N);
        end
        function ComputeBoundaryValuesSubRef(obj, DiscretizationModel, Formulation, CellsSelected, CellsSelected_old)
            % merge the new accepted fractional flow with those occepted at
            % the previous sub-ref
            
            CellsSelected.f = Formulation.f .* (CellsSelected_old.ActCells) + ...
                CellsSelected_old.f .* (1 - CellsSelected_old.ActCells);

            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N  = DiscretizationModel.ReservoirGrid.N;
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
        end
        
        function SetCorrectActiveCells(obj, CellsSelected)
            obj.f = CellsSelected.f;
            obj.ActCells = CellsSelected.ActCells;
            obj.BCFluxMatrix = CellsSelected.BCFluxMatrix;
        end
        function ViscousMatrixLTS(obj, Grid, Formulation)
            %Builds Upwind Flux matrix
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;                                   
            q = min(Formulation.Qwells, 0) .* obj.ActCells;  
             
            Utot = Formulation.Utot;
            % right to left and top to bottom (negative x, y, z)
            Xneg = min(Utot.x, 0); 
            Yneg = min(Utot.y, 0);
            Zneg = min(Utot.z, 0);
            % make them vectors 
            x1 = reshape(Xneg(1:Nx,:,:),N,1);
            y1 = reshape(Yneg(:,1:Ny,:),N,1);
            z1 = reshape(Zneg(:,:,1:Nz),N,1);
            
            ActFluxes = Formulation.MatrixAssembler.ActInterfaces;
            
            %select only the active fluxes
            x1A = reshape(ActFluxes.x(1:Nx,:,:),N,1);
            y1A = reshape(ActFluxes.y(:,1:Ny,:),N,1);
            z1A = reshape(ActFluxes.z(:,:,1:Nz),N,1);
            
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
            x2A = reshape(ActFluxes.x(2:Nx+1,:,:),N,1);
            y2A = reshape(ActFluxes.y(:,2:Ny+1,:),N,1);
            z2A = reshape(ActFluxes.z(:,:,2:Nz+1),N,1);
            
            x2 = x2 .* x2A;
            y2 = y2 .* y2A;
            z2 = z2 .* z2A;
            
            % Assemble matrix
            DiagVecs = [z2, y2, x2, +x1-x2+y1-y2+z1-z2+q, -x1, -y1, -z1]; % diagonal vectors
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny]; % diagonal index
            obj.Vr = spdiags(DiagVecs, DiagIndx, N, N);
        end
        function Residual = AddBC2Residual(obj, Residual, ProductionSystem, Formulation, DiscretizationModel, State0, dt)
            
            obj.ViscousMatrixLTS(DiscretizationModel.ReservoirGrid, Formulation);   
            
            rho = ProductionSystem.Reservoir.State.Properties('rho_1').Value;
            % (1) Put to zero the accepted residual cells 
            Residual = (Residual + Formulation.V * (Formulation.f .* rho)) .* obj.ActCells; 
            % (2) Add the contribution of the viscous matrix 
            Residual = Residual - obj.Vr * (Formulation.f .* rho);
            % (3) Add the contribution of the accelpted fluxes at the
            % boundary
            Residual = Residual - (obj.BCFluxMatrix * (obj.f .* rho)).* obj.ActCells ;

        end
        function Jacobian = AddBC2Jacobian(obj, Jacobian, ProductionSystem, Formulation, DiscretizationModel, dt)
            N = DiscretizationModel.ReservoirGrid.N;
            rho = ProductionSystem.Reservoir.State.Properties('rho_1').Value;
            
            Jacobian = (Jacobian + Formulation.V * spdiags(Formulation.df.*rho,0,N,N)).*spdiags(obj.ActCells,0,N,N) + ...
                - obj.Vr * spdiags(Formulation.df .* rho .* obj.ActCells,0,N,N); 
        end
    end
end
