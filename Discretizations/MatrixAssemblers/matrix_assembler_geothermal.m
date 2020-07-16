% Matrix assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef matrix_assembler_geothermal < matrix_assembler
    properties
        
    end
    methods
        function [Tph, Gph, Thph, Tk] = TransmissibilityMatrix(obj, Grid, Upwind, Mob, rho, h, RhoInt)
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
                    
                    %% Heat Convection Flux Transmissibility Matrix
                    % Apply Upwind Operator
                    Mob_rho_h_Upwind = Upwind * (Mob.*rho.*h);
                    % Transmisibility Matrix
                    Thph = Grid.ConnectivityMatrix * spdiags(Grid.Trans.*Mob_rho_h_Upwind,0,N_Faces,N_Faces) * Grid.ConnectivityMatrix';
                    
                    %% Heat Conduction Flux Transmissibility matrix
                    Tk = Grid.ConnectivityMatrix * spdiags(Grid.HeatTrans                ,0,N_Faces,N_Faces) * Grid.ConnectivityMatrix';
                    
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
                    Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
                    Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
                    Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
                    Tph = obj.ReshapeCartesianTransmisibility(Nx, Ny, Nz, N, Tx, Ty, Tz);
                    
                    % Gravity Matrix
                    Tx(2:Grid.Nx,:,:)= Tx(2:Grid.Nx,:,:) .* RhoInt.x(2:Grid.Nx,:,:);
                    Ty(:,2:Grid.Ny,:)= Ty(:,2:Grid.Ny,:) .* RhoInt.y(:,2:Grid.Ny,:);
                    Tz(:,:,2:Grid.Nz)= Tz(:,:,2:Grid.Nz) .* RhoInt.z(:,:,2:Grid.Nz);
                    Gph = obj.ReshapeCartesianTransmisibility(Nx, Ny, Nz, N, Tx, Ty, Tz);
                    
                    %% Heat Convection Flux Transmissibility Matrix
                    Tx = zeros(Nx+1, Ny, Nz);
                    Ty = zeros(Nx, Ny+1, Nz);
                    Tz = zeros(Nx, Ny, Nz+1);
                    
                    % Apply Upwind Operator
                    Mupx = Upwind.x*(Mob .* rho .* h);
                    Mupy = Upwind.y*(Mob .* rho .* h);
                    Mupz = Upwind.z*(Mob .* rho .* h);
                    Mupx = reshape(Mupx, Nx, Ny, Nz);
                    Mupy = reshape(Mupy, Nx, Ny, Nz);
                    Mupz = reshape(Mupz, Nx, Ny, Nz);
                    
                    Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
                    Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
                    Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
                    Thph = obj.ReshapeCartesianTransmisibility(Nx, Ny, Nz, N, Tx, Ty, Tz);
                    
                    %% Heat Conduction Flux Transmissibility matrix
                    THx = zeros(Nx+1, Ny, Nz);
                    THy = zeros(Nx, Ny+1, Nz);
                    THz = zeros(Nx, Ny, Nz+1);
                    
                    THx(2:Nx,:,:)= Grid.THx(2:Nx,:,:);
                    THy(:,2:Ny,:)= Grid.THy(:,2:Ny,:);
                    THz(:,:,2:Nz)= Grid.THz(:,:,2:Nz);
                    Tk = obj.ReshapeCartesianTransmisibility(Nx, Ny, Nz, N, THx, THy, THz); % Transmisibility of rock conductivity
            end
        end
    end
end
