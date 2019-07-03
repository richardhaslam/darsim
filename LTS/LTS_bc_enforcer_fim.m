% LTS bc enforcer for FIM coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_bc_enforcer_fim < LTS_bc_enforcer
    properties
       ActCells
       BCFluxMatrix
       DirBCFlux
       CellBc
       BcNeighbour
       T
    end
    methods
        function ComputeBoundaryFluxes(obj, DiscretizationModel, Formulation, ActCells)
           
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            Nz = DiscretizationModel.ReservoirGrid.Nz;
            N  = DiscretizationModel.ReservoirGrid.N;
            
            BC.x = zeros(Nx+1,Ny,Nz);
            BC.y = zeros(Nx,Ny+1,Nz);
            BC.z = zeros(Nx,Ny,Nz+1);
            obj.ActCells = ActCells;
            ActCellsM = reshape(obj.ActCells, Nx, Ny, Nz);      
            
            U = Formulation.Utot;
            
            %fluxes at the boundaty between an accepted and a refuced cell
            BC.x(2:Nx,:,:) = (ActCellsM(1:Nx-1,:,:) + ActCellsM(2:Nx,:,:) == 1);
            BC.y(:,2:Ny,:) = (ActCellsM(:,1:Ny-1,:) + ActCellsM(:,2:Ny,:) == 1);
            BC.z(:,:,2:Nz) = (ActCellsM(:,:,1:Nz-1) + ActCellsM(:,:,2:Nz) == 1);
            
            %fluxes at the boundaty between an accepted and a refuced cell
            BCFluxes.x = U.x .* BC.x;
            BCFluxes.y = U.y .* BC.y;
            BCFluxes.z = U.z .* BC.z;
            
            % right to left and top to bottom (negative x, y, z)
            Xneg = min(BCFluxes.x, 0);
            Yneg = min(BCFluxes.y, 0);
            Zneg = min(BCFluxes.z, 0);
            
            x1 = reshape(Xneg(1:Nx,:,:),N,1);
            y1 = reshape(Yneg(:,1:Ny,:),N,1);
            z1 = reshape(Zneg(:,:,1:Nz),N,1);
            
            % left to right and bottom to top (positive x, y, z)
            Xpos = max(BCFluxes.x, 0);
            Ypos = max(BCFluxes.y, 0);
            Zpos = max(BCFluxes.z, 0);
            
            x2 = reshape(Xpos(2:Nx+1,:,:), N, 1);
            y2 = reshape(Ypos(:,2:Ny+1,:), N, 1);
            z2 = reshape(Zpos(:,:,2:Nz+1), N, 1);
            
            DiagVecs = [z2, y2, x2, +x1-x2+y1-y2+z1-z2, -x1, -y1, -z1]; % diagonal vectors
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            
            obj.BCFluxMatrix = spdiags(DiagVecs, DiagIndx, N, N);
            
            % Let us choose one cell in x to have a Dirichlet b.c.
            [i, j, k] = ind2sub(size(BC.x),find(BC.x));
            if i(end) < Nx
                Index = i(end) - 1 + (j(end)-1) * Nx + (k(end) -1) * Nx * Ny;
                neighbour = Index + 1;
            else
                Index = i(1) - 1 + (j(1)-1) * Nx + (k(1) -1) * Nx * Ny;
                neighbour = Index + 1;
            end
            cellbc = Index;
            obj.CellBc  = cellbc;
            obj.BcNeighbour = neighbour;
            obj.T = DiscretizationModel.ReservoirGrid.dy * DiscretizationModel.ReservoirGrid.dz / DiscretizationModel.ReservoirGrid.dx;
        end
        function Residual = AddBC2Residual(obj, Residual, ProductionSystem, Formulation)
            N = length(Formulation.Mob(:,1));
            f = zeros(N, 2);
            f(:, 1) = Formulation.Mob(:,1) ./ (Formulation.Mob(:,1) + Formulation.Mob(:, 2));
            f(:, 2) = 1 - f(:, 1); 
            
            Dirichlet = ones(N,1);
            Dirichlet(obj.CellBc) = 0;             
            KH = ProductionSystem.Reservoir.K(obj.CellBc,1);
             
            End = 0;
            for ph=1:2
                %% Add flux b.c.
                rho = ProductionSystem.Reservoir.State.Properties(strcat('rho_', num2str(ph)));
                Start = End + 1;
                End = Start -  1 + N;
                Residual(Start:End) = Residual(Start:End).* obj.ActCells - (obj.BCFluxMatrix * f(:, ph) .* rho.Value) .* obj.ActCells .* Dirichlet;
                
                %% Add 1 pressure b.c.
                P = ProductionSystem.Reservoir.State.Properties(strcat('P_', num2str(ph)));
                Index = obj.CellBc + (ph-1)*N;
                Residual(Index) = Residual(Index) - obj.T * KH * Formulation.Mob(obj.CellBc, ph) * rho.Value(obj.CellBc) * (P.Value(obj.BcNeighbour)  - P.Value(obj.CellBc)); 
            end
        end
        function Jacobian = AddBC2Jacobian(obj, Jacobian, ProductionSystem, Formulation, DiscretizationModel)
            % Add pressure b.c.
            N = DiscretizationModel.ReservoirGrid.N;
            Mob = Formulation.Mob;
            dMob = Formulation.dMob;
            rho1 = ProductionSystem.Reservoir.State.Properties('rho_1');  
            rho2 = ProductionSystem.Reservoir.State.Properties('rho_2');  
            KH = ProductionSystem.Reservoir.K(obj.CellBc,1);
            P1 = ProductionSystem.Reservoir.State.Properties('P_1');  
            P2 = ProductionSystem.Reservoir.State.Properties('P_2');  
            
            % row and column indexes
            I1 = obj.CellBc;
            I2 = obj.CellBc + N;
            Ip = obj.CellBc;
            IS = obj.CellBc + N;
            % fill in the matrix
            Jacobian(I1, Ip) = Jacobian(I1, Ip) + obj.T * KH *...
                Mob(obj.CellBc, 1) * rho1.Value(obj.CellBc);
            Jacobian(I1, IS) = Jacobian(I1, IS) - obj.T * KH *...
                dMob(obj.CellBc, 1) * rho1.Value(obj.CellBc) * (P1.Value(obj.BcNeighbour)  - P1.Value(obj.CellBc)) ; 
            Jacobian(I2, Ip) = Jacobian(I2, Ip) + obj.T * KH *...
                Mob(obj.CellBc, 2) * rho2.Value(obj.CellBc);
            Jacobian(I2, IS) = Jacobian(I2, IS) - obj.T * KH *...
                dMob(obj.CellBc, 2) * rho2.Value(obj.CellBc) * (P2.Value(obj.BcNeighbour)  - P2.Value(obj.CellBc)) ;
        end
    end
end
