% Matrix assembler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef matrix_assembler_geothermal < matrix_assembler
    % This is a subclass of the matrix_assembler class, therefore you only
    % have to define the HeatTransmissibilityMatrix() function as the
    % regular one is in the super-class
    
    properties
        
    end
    methods
        function [Th, Gh] = ConvectiveHeatTransmissibilityMatrix(obj, Grid, UpWind, Mob, rho, h, RhoInt)
            % Th : Heat convection transmissibility (Tph*rho)
            % Gh : Heat convection gravity
            
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            
            %% Th Matrix
            Tx = zeros(Nx+1, Ny, Nz);
            Ty = zeros(Nx, Ny+1, Nz);
            Tz = zeros(Nx, Ny, Nz+1);

            % Apply upwind operator; Phase enthalpy is included in Upwind
            Mupx = UpWind.x*(Mob .* rho .* h);
            Mupy = UpWind.y*(Mob .* rho .* h);
            Mupz = UpWind.z*(Mob .* rho .* h);
            Mupx = reshape(Mupx, Nx, Ny, Nz);
            Mupy = reshape(Mupy, Nx, Ny, Nz);
            Mupz = reshape(Mupz, Nx, Ny, Nz);

            % Transmissibility matrix
            Tx(2:Nx,:,:)= Grid.Tx(2:Nx,:,:).*Mupx(1:Nx-1,:,:);
            Ty(:,2:Ny,:)= Grid.Ty(:,2:Ny,:).*Mupy(:,1:Ny-1,:);
            Tz(:,:,2:Nz)= Grid.Tz(:,:,2:Nz).*Mupz(:,:,1:Nz-1);
            Th = obj.ReshapeTransmisibility(Grid, Tx, Ty, Tz); 
            
            % Gravity Matrix
            Tx(2:Grid.Nx,:,:)= Tx(2:Grid.Nx,:,:) .* RhoInt.x(2:Grid.Nx,:,:);
            Ty(:,2:Grid.Ny,:)= Ty(:,2:Grid.Ny,:) .* RhoInt.y(:,2:Grid.Ny,:);
            Tz(:,:,2:Grid.Nz)= Tz(:,:,2:Grid.Nz) .* RhoInt.z(:,:,2:Grid.Nz);
            Gh = obj.ReshapeTransmissibility(Grid, Tx, Ty, Tz);

        end
        function Tk = ConductiveHeatTransmissibilityMatrix(obj, Grid)
            % Tk : Heat conduction transmissibility
            
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;

            %% Tk Matrix
            THx = zeros(Nx+1, Ny, Nz);
            THy = zeros(Nx, Ny+1, Nz);
            THz = zeros(Nx, Ny, Nz+1);
            
            THx(2:Nx,:,:)= Grid.THx(2:Nx,:,:); 
            THy(:,2:Ny,:)= Grid.THy(:,2:Ny,:);
            THz(:,:,2:Nz)= Grid.THz(:,:,2:Nz);
            Tk = obj.ReshapeTransmisibility(Grid, THx, THy, THz); % Transmisibility of rock conductivity
        end
    end
end

