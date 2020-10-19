% Phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 26 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef phase < matlab.mixin.Heterogeneous & handle
    properties
        mu % Reference Viscosity
        sr % Irriducible saturation
    end
    methods
        function [Upwind, U] = UpWindAndRockFluxes(obj, Grid, P, RhoInt)
            switch class(Grid)
                case('corner_point_grid')
                    Neighbor1Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index;
                    Neighbor2Index = Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index;
                    
                    % Gravitational velocities
                    Ug = - ( Grid.Depth(Neighbor2Index) - Grid.Depth(Neighbor1Index) ) .* RhoInt;
                    
                    % Compute 'rock' fluxes ([m^3/s])
                    U = - Grid.Trans .* ( P(Neighbor2Index)-P(Neighbor1Index) - Ug );
                    
                    % Use velocity to build upwind operator
                    nc = Grid.N;
                    nf = length(Grid.Trans);
                    C = (U>=0) .* Grid.CornerPointGridData.Internal_Faces.CellNeighbor1Index + ...
                        (U< 0) .* Grid.CornerPointGridData.Internal_Faces.CellNeighbor2Index;
                    Upwind = sparse( (1:nf)', C , ones(nf,1) , nf , nc );
                    
                case('cartesian_grid')
                    Nx = Grid.Nx;
                    Ny = Grid.Ny;
                    Nz = Grid.Nz;
                    N = Grid.N;

                    % Gravitational velocities
                    depth = reshape(Grid.Depth, Nx, Ny, Nz);
                    Ugx = zeros(Nx+1, Ny, Nz);
                    Ugy = zeros(Nx, Ny+1, Nz);
                    Ugz = zeros(Nx, Ny, Nz+1);
                    Ugx(2:Nx,:,:) = (depth(1:Nx-1,:,:) - depth(2:Nx,:,:)) .* RhoInt.x(2:Nx,:,:);
                    Ugy(:,2:Ny,:) = (depth(:,1:Ny-1,:) - depth(:,2:Ny,:)) .* RhoInt.y(:,2:Ny,:);
                    Ugz(:,:,2:Nz) = (depth(:,:,1:Nz-1) - depth(:,:,2:Nz)) .* RhoInt.z(:,:,2:Nz);
                    
                    P = reshape(P, Nx, Ny, Nz);
                    %% Compute 'rock' fluxes ([m^3/s])
                    U.x = zeros(Nx+1,Ny,Nz);
                    U.y = zeros(Nx,Ny+1,Nz);
                    U.z = zeros(Nx,Ny,Nz+1);
                    U.x(2:Nx,:,:) = Grid.Tx(2:Nx,:,:) .* ( (P(1:Nx-1,:,:)-P(2:Nx,:,:)) - Ugx(2:Nx,:,:) );
                    U.y(:,2:Ny,:) = Grid.Ty(:,2:Ny,:) .* ( (P(:,1:Ny-1,:)-P(:,2:Ny,:)) - Ugy(:,2:Ny,:) );
                    U.z(:,:,2:Nz) = Grid.Tz(:,:,2:Nz) .* ( (P(:,:,1:Nz-1)-P(:,:,2:Nz)) - Ugz(:,:,2:Nz) );
                    
                    %% Use velocity to build upwind operator
                    L = reshape((U.x(2:Nx+1,:,:) >= 0), N, 1);
                    R = reshape((U.x(1:Nx,:,:) < 0), N, 1);
                    B = reshape((U.y(:,2:Ny+1,:) >= 0), N, 1);
                    T = reshape((U.y(:,1:Ny,:) < 0), N, 1);
                    Down = reshape((U.z(:,:,2:Nz+1) >= 0), N, 1);
                    Up = reshape((U.z(:,:,1:Nz) < 0), N, 1);
                    
                    DiagVecs = [L, R];
                    DiagIndx = [0, 1];
                    Upwind.x = spdiags(DiagVecs, DiagIndx, N, N);
                    DiagVecs = [B, T];
                    DiagIndx = [0, Nx];
                    Upwind.y = spdiags(DiagVecs, DiagIndx, N, N);
                    DiagVecs = [Down, Up];
                    DiagIndx = [0, Nx*Ny];
                    Upwind.z = spdiags(DiagVecs, DiagIndx, N, N);
            end
        end
    end
end