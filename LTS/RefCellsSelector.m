%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef RefCellsSelector < handle
    properties
        tol = 1e-2
        NSten = 1
        ActFluxes
        ActCells
        ViscousMatrixValue
        f
        df
    end
    methods
        function SelectRefCells(obj, ProductionSystem, Grid, Formulation)
            
            Utot = Formulation.Utot;
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            N = Grid.N;
            % We choose to resolve with a fine dt all cells which were
            % flooded during the coarse dt. 
            %For now it only works for the reservoir and not for the
            %fractures
            dS = ProductionSystem.Reservoir.State.Properties('S_1').Value - ...
                 ProductionSystem.Reservoir.State_old.Properties('S_1').Value;
            % vector of active cells
            obj.ActCells = abs(dS) >= obj.tol;
            
            % to add also the cell newt to the upwind velocity
            
            ActCellsM = reshape(obj.ActCells, Nx, Ny, Nz);
            
            % if you want to increase the set just apply a for cycle
            for i = 1: obj.NSten
                ActCellsM(2:Nx,:,:) = ActCellsM(2:Nx,:,:) + Utot.x(2:Nx,:,:) .* ActCellsM(1:Nx-1,:,:);
                ActCellsM(1:Nx-1,:,:) = ActCellsM(1:Nx-1,:,:) + Utot.x(1:Nx-1,:,:) .*  ActCellsM(2:Nx,:,:);
            
                ActCellsM(:,2:Ny,:) = ActCellsM(:,2:Ny,:) + Utot.y(:,2:Ny,:) .* ActCellsM(:,1:Ny-1,:);
                ActCellsM(:,1:Ny-1,:) = ActCellsM(:,1:Ny-1,:) + Utot.y(:,1:Ny-1,:) .* ActCellsM(:,2:Ny,:);
            
                ActCellsM(:,:,2:Nz) = ActCellsM(:,:,2:Nz) + Utot.z(:,:,2:Nz) .* ActCellsM(:,:,1:Nz-1);
                ActCellsM(:,:,1:Nz-1) = ActCellsM(:,:,1:Nz-1) + Utot.z(:,:,1:Nz-1) .* ActCellsM(:,:,2:Nz);
            end
%             ActCellsM(2:Nx,:,:) = ActCellsM(2:Nx,:,:) + ActCellsM(1:Nx-1,:,:); 
%             ActCellsM(:,2:Ny,:) = ActCellsM(:,2:Ny,:) + ActCellsM(:,1:Ny-1,:);
            obj.ActCells = reshape(ActCellsM >0, N, 1);
            
        end
        function ActFluxes = SelectRefFluxes(obj, Grid)
            
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            % From the ActCells mask we reconstruct the good fluxes.
            ActFluxes.x = zeros(Nx+1,Ny,Nz);
            ActFluxes.y = zeros(Nx,Ny+1,Nz);
            ActFluxes.z = zeros(Nx,Ny,Nz+1);
               
            ActCellsM = reshape(obj.ActCells, Nx, Ny, Nz);            
            
            ActFluxes.x(1,:,:)    = ActCellsM(1,:,:);
            ActFluxes.x(Nx+1,:,:) = ActCellsM(Nx,:,:);
            ActFluxes.x(2:Nx,:,:) = (ActCellsM(1:Nx-1,:,:) + ActCellsM(2:Nx,:,:)) == 2;
            
            ActFluxes.y(:,1,:)    = ActCellsM(:,1,:);
            ActFluxes.y(:,Ny+1,:) = ActCellsM(:,Ny,:);
            ActFluxes.y(:,2:Ny,:) = (ActCellsM(:,1:Ny-1,:) + ActCellsM(:,2:Ny,:)) == 2;
            
            ActFluxes.z(:,:,1)    = ActCellsM(:,:,1);
            ActFluxes.z(:,:,Nz+1) = ActCellsM(:,:,Nz);
            ActFluxes.z(:,:,2:Nz) = (ActCellsM(:,:,1:Nz-1) + ActCellsM(:,:,2:Nz)) == 2; 
            obj.ActFluxes = ActFluxes;
        end
        function ComputeBoundaryValues(obj, DiscretizationModel, Formulation)
            % store the past values of the fractional flow and the
            % derivatives of the fractional flow
            obj.f = Formulation.f;
            obj.df = Formulation.df;
            
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N  = DiscretizationModel.ReservoirGrid.N;
            
            BCFluxes.x = zeros(Nx+1,Ny,Nz);
            BCFluxes.y = zeros(Nx,Ny+1,Nz);
            BCFluxes.z = zeros(Nx,Ny,Nz+1);
            ActCellsM = reshape(obj.ActCells, Nx, Ny, Nz);
            
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

            obj.ViscousMatrixValue = spdiags(DiagVecs, DiagIndx, N, N);
        end
    end
end