% Compositional Formulation base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 October 2016
%Last modified: 7 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Compositional_formulation < formulation
    properties
        NofComponents
    end
    methods
        function obj = Compositional_formulation(n_components)
            obj@formulation();
            obj.NofComponents = n_components;
            obj.Tph = cell(n_components, 2);
            obj.Gph = cell(n_components, 2);
        end
        function  TransmissibilityMatrix(obj, Grid, Rho, RhoInt, x, i)
            %%%Transmissibility matrix construction
            Tx = zeros(Grid.Nx+1, Grid.Ny, Grid.Nz);
            Ty = zeros(Grid.Nx, Grid.Ny+1, Grid.Nz);
            Tz = zeros(Grid.Nx, Grid.Ny, Grid.Nz+1);
            
            %% PHASE 1 
            %Apply upwind operator
            Mupx = obj.UpWind(1).x*(obj.Mob(:,1) .* Rho(:,1) .* x(:,1)); 
            Mupy = obj.UpWind(1).y*(obj.Mob(:,1) .* Rho(:,1) .* x(:,1));
            Mupz = obj.UpWind(1).z*(obj.Mob(:,1) .* Rho(:,1) .* x(:,1));
            Mupx = reshape(Mupx, Grid.Nx, Grid.Ny, Grid.Nz);
            Mupy = reshape(Mupy, Grid.Nx, Grid.Ny, Grid.Nz);
            Mupz = reshape(Mupz, Grid.Nx, Grid.Ny, Grid.Nz);
            
            Tx(2:Grid.Nx,:,:)= Grid.Tx(2:Grid.Nx,:,:).*Mupx(1:Grid.Nx-1,:,:);
            Ty(:,2:Grid.Ny,:)= Grid.Ty(:,2:Grid.Ny,:).*Mupy(:,1:Grid.Ny-1,:);
            Tz(:,:,2:Grid.Nz)= Grid.Tz(:,:,2:Grid.Nz).*Mupz(:,:,1:Grid.Nz-1);
            
            %Construct matrix
            x1 = reshape(Tx(1:Grid.Nx,:,:), Grid.N, 1);
            x2 = reshape(Tx(2:Grid.Nx+1,:,:), Grid.N, 1);
            y1 = reshape(Ty(:,1:Grid.Ny,:), Grid.N, 1);
            y2 = reshape(Ty(:,2:Grid.Ny+1,:), Grid.N, 1);
            z1 = reshape(Tz(:,:,1:Grid.Nz), Grid.N, 1);
            z2 = reshape(Tz(:,:,2:Grid.Nz+1), Grid.N, 1);
            
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
            obj.Tph{i,1} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
            
            % Gravity Matrix
            Tx(2:Grid.Nx,:,:)= Tx(2:Grid.Nx,:,:) .* RhoInt(1).x(2:Grid.Nx,:,:);
            Ty(:,2:Grid.Ny,:)= Ty(:,2:Grid.Ny,:) .* RhoInt(1).y(:,2:Grid.Ny,:);
            Tz(:,:,2:Grid.Nz)= Tz(:,:,2:Grid.Nz) .* RhoInt(1).z(:,:,2:Grid.Nz);
            
            %Construct matrix
            x1 = reshape(Tx(1:Grid.Nx,:,:), Grid.N, 1);
            x2 = reshape(Tx(2:Grid.Nx+1,:,:), Grid.N, 1);
            y1 = reshape(Ty(:,1:Grid.Ny,:), Grid.N, 1);
            y2 = reshape(Ty(:,2:Grid.Ny+1,:), Grid.N, 1);
            z1 = reshape(Tz(:,:,1:Grid.Nz), Grid.N, 1);
            z2 = reshape(Tz(:,:,2:Grid.Nz+1), Grid.N, 1);
            
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
            obj.Gph{i,1} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
            
            %% PHASE 2 
            %Apply upwind operator
            Mupx = obj.UpWind(2).x*(obj.Mob(:,2) .* Rho(:,2) .* x(:,2));
            Mupy = obj.UpWind(2).y*(obj.Mob(:,2) .* Rho(:,2) .* x(:,2));
            Mupz = obj.UpWind(2).z*(obj.Mob(:,2) .* Rho(:,2) .* x(:,2));
            Mupx = reshape(Mupx, Grid.Nx, Grid.Ny, Grid.Nz);
            Mupy = reshape(Mupy, Grid.Nx, Grid.Ny, Grid.Nz);
            Mupz = reshape(Mupz, Grid.Nx, Grid.Ny, Grid.Nz);
            
            Tx(2:Grid.Nx,:,:) = Grid.Tx(2:Grid.Nx,:,:).*Mupx(1:Grid.Nx-1,:,:);
            Ty(:,2:Grid.Ny,:) = Grid.Ty(:,2:Grid.Ny,:).*Mupy(:,1:Grid.Ny-1,:);
            Tz(:,:,2:Grid.Nz) = Grid.Tz(:,:,2:Grid.Nz).*Mupz(:,:,1:Grid.Nz-1);
            
            %Construct matrix
            x1 = reshape(Tx(1:Grid.Nx,:,:), Grid.N, 1);
            x2 = reshape(Tx(2:Grid.Nx+1,:,:), Grid.N, 1);
            y1 = reshape(Ty(:,1:Grid.Ny,:), Grid.N, 1);
            y2 = reshape(Ty(:,2:Grid.Ny+1,:), Grid.N, 1);
            z1 = reshape(Tz(:,:,1:Grid.Nz), Grid.N, 1);
            z2 = reshape(Tz(:,:,2:Grid.Nz+1), Grid.N, 1);
            
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
            obj.Tph{i,2} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
            
            
            % Gravity Matrix
            Tx(2:Grid.Nx,:,:) = Tx(2:Grid.Nx,:,:) .* RhoInt(2).x(2:Grid.Nx,:,:);
            Ty(:,2:Grid.Ny,:) = Ty(:,2:Grid.Ny,:) .* RhoInt(2).y(:,2:Grid.Ny,:);
            Tz(:,:,2:Grid.Nz) = Tz(:,:,2:Grid.Nz) .* RhoInt(2).z(:,:,2:Grid.Nz);
            
            %Construct matrix
            x1 = reshape(Tx(1:Grid.Nx,:,:), Grid.N, 1);
            x2 = reshape(Tx(2:Grid.Nx+1,:,:), Grid.N, 1);
            y1 = reshape(Ty(:,1:Grid.Ny,:), Grid.N, 1);
            y2 = reshape(Ty(:,2:Grid.Ny+1,:), Grid.N, 1);
            z1 = reshape(Tz(:,:,1:Grid.Nz), Grid.N, 1);
            z2 = reshape(Tz(:,:,2:Grid.Nz+1), Grid.N, 1);
            
            DiagVecs = [-z2,-y2,-x2,z2+y2+x2+y1+x1+z1,-x1,-y1,-z1];
            DiagIndx = [-Grid.Nx*Grid.Ny, -Grid.Nx, -1, 0, 1, Grid.Nx, Grid.Nx*Grid.Ny];
            obj.Gph{i,2} = spdiags(DiagVecs, DiagIndx, Grid.N, Grid.N);
            
        end
        function q = ComputeSourceTerms(obj, N, Wells)
            q = zeros(N, obj.NofComponents);
            %Injectors
            for i=1:Wells.NofInj
                c = Wells.Inj(i).Cells;
                q(c, :) = Wells.Inj(i).QComponents(:,:);
            end
            %Producers
            for i=1:Wells.NofProd
                c = Wells.Prod(i).Cells;
                q(c, :) = Wells.Prod(i).QComponents(:,:);
            end
        end
        function AverageMassOnCoarseBlocks(obj, Status, FluidModel, R, P)
            % Perform Average for ADM
            z_rest = R * Status.z;
            Status.z = P * z_rest;
            % Update other unknwons as well 
            obj.UpdatePhaseCompositions(Status, FluidModel);
        end
        function UpdatePhaseCompositions(obj, Status, FluidModel)
            %% 2. Perform composition update ((x, ni) = f(T, p, z))
            obj.SinglePhase = FluidModel.Flash(Status);
            
            %% 3. Compute Densities (rho = rho(p, x, T))
            FluidModel.ComputePhaseDensities(Status);
            
            %% 4. Compute Saturations (S = S(z, x))
            FluidModel.ComputePhaseSaturation(Status, obj.SinglePhase);
            
            %% 5. Compute Total Density (rhoT = rhoT(S, rho))
            FluidModel.ComputeTotalDensity(Status);
            
            %% 6. Compute Pc (Pc = Pc(S))
            FluidModel.ComputePc(Status);
        end
    end
end