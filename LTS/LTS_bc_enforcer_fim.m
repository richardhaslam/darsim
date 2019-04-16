%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_bc_enforcer_fim < LTS_bc_enforcer
    properties
                
    end
    methods
        function ComputeBoundaryValues(obj, DiscretizationModel, Formulation, ActCells)
           
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N  = DiscretizationModel.ReservoirGrid.N;
            
            BCFluxes.x = zeros(Nx+1,Ny,Nz);
            BCFluxes.y = zeros(Nx,Ny+1,Nz);
            BCFluxes.z = zeros(Nx,Ny,Nz+1);
            ActCellsM = reshape(ActCells, Nx, Ny, Nz);
            
            for ph = 1:2
                % Compute the fluxes values
                U = Formulation.U{ph, 1};
                
                % right to left and top to bottom (negative x, y, z)
                Xneg = min(Utot.x, 0);
                Yneg = min(Utot.y, 0);
                Zneg = min(Utot.z, 0);
                
                x1 = reshape(Xneg(1:Nx,:,:),N,1);
                y1 = reshape(Yneg(:,1:Ny,:),N,1);
                z1 = reshape(Zneg(:,:,1:Nz),N,1);
                
                % left to right and bottom to top (positive x, y, z)
                Xpos = max(U.x, 0);
                Ypos = max(U.y, 0);
                Zpos = max(U.z, 0);
                
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
                
                obj.ViscousMatrixValue{ph, 1} = spdiags(DiagVecs, DiagIndx, N, N);
            end
        end
        function ComputeBoundaryValuesSubRef(obj, DiscretizationModel, Formulation, CellsSelected_old, ActCells)
            % merge the new accepted fractional flow with those occepted at
            % the previous sub-ref
            
            obj.f = Formulation.f .* (CellsSelected_old.ActCells) + ...
                CellsSelected_old.BCEnforcer.f .* (1 - CellsSelected_old.ActCells);
            
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N  = DiscretizationModel.ReservoirGrid.N;
            % merge the new fluxes interfaces with those obtained at the
            % previous sub-ref.
            
            BCFluxes.x = zeros(Nx+1,Ny,Nz);
            BCFluxes.y = zeros(Nx,Ny+1,Nz);
            BCFluxes.z = zeros(Nx,Ny,Nz+1);
            ActCellsM = reshape(ActCells, Nx, Ny, Nz);
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
            
            %fluxes at the boundary between an accepted and a refuced cell
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
            
            obj.ViscousMatrixValue = spdiags(DiagVecs, DiagIndx, N, N);
        end    
    end
end
