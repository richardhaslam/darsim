%  Basis functions updater
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 16 August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef bf_updater_bilin < bf_updater
    properties
        
    end
    methods 
         function ConstructPressureSystem(obj, FineGrid, K, s, FluidModel)
            % Initialize local variables
            Nx = FineGrid.Nx;
            Ny = FineGrid.Ny;
            N = FineGrid.N;
            dx = FineGrid.dx;
            dy = FineGrid.dy;
            Ax = FineGrid.Ax;
            Ay = FineGrid.Ay;
            Mob = FluidModel.ComputePhaseMobilities(s);
            Mobt = sum(Mob,2);
            
            K(:, 1) = mean(K(:,1)) .* Mobt;
            K(:, 2) = mean(K(:,2)) .* Mobt;
            
            % Harmonic average of permeability
            kx = reshape(K(:,1), Nx, Ny);
            ky = reshape(K(:,2), Nx, Ny);
            Kx = zeros(Nx+1, Ny);
            Ky = zeros(Nx, Ny+1);
            Kx(2:Nx,:) = 2*kx(1:Nx-1,:) .* kx(2:Nx,:) ./ (kx(1:Nx-1,:) + kx(2:Nx,:));
            Ky(:,2:Ny) = 2*ky(:,1:Ny-1) .* ky(:,2:Ny) ./ (ky(:,1:Ny-1) + ky(:,2:Ny));
            
            % Transmissibility
            Tx = zeros(Nx+1, Ny);
            Ty = zeros(Nx, Ny+1);
            Tx(2:Nx,:) = Ax./dx .* Kx(2:Nx,:);
            Ty(:,2:Ny) = Ay./dy .* Ky(:,2:Ny);
            
            %Construct pressure matrix
            x1 = reshape(Tx(1:Nx,:),N,1);
            x2 = reshape(Tx(2:Nx+1,:),N,1);
            y1 = reshape(Ty(:,1:Ny),N,1);
            y2 = reshape(Ty(:,2:Ny+1),N,1);
            DiagVecs = [-y2,-x2,y2+x2+y1+x1,-x1,-y1];
            DiagIndx = [-Nx,-1,0,1,Nx];
            obj.A = spdiags(DiagVecs,DiagIndx,N,N);
        end
    end
end