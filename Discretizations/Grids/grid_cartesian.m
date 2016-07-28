%  Cartesian grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 14 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef grid_cartesian < grid
    properties
        Nx
        Ny
        Nz
        dx
        dy
        dz
        Ax
        Ay
        Volume
        Tx
        Ty
        I
        J
    end
    methods
        function obj = grid_cartesian(nx, ny, nz)
            obj.Nx = nx;
            obj.Ny = ny;
            obj.Nz = nz;
            obj.N = nx * ny * nz;
            obj.Active = ones(obj.N, 1);
            obj.Tx = zeros(nx+1, ny);
            obj.Ty = zeros(nx, ny+1);
        end
        function Initialize(obj, Reservoir)
            obj.dx = Reservoir.Length/obj.Nx;
            obj.dy = Reservoir.Width/obj.Ny;
            obj.dz = Reservoir.Depth/obj.Nz;
            obj.Ax = obj.dy * obj.dz;
            obj.Ay = obj.dx * obj.dz;
            obj.Volume = obj.dx * obj.dy * obj.dz;
            obj.ComputeRockTransmissibilities(Reservoir.K);
        end
        function ComputeRockTransmissibilities(obj, K)
            %Harmonic average of permeability.
            Kx = reshape(K(:,1), obj.Nx, obj.Ny);
            Ky = reshape(K(:,2), obj.Nx, obj.Ny);
            
            KHx = zeros(obj.Nx+1, obj.Ny);
            KHy = zeros(obj.Nx, obj.Ny+1);
            KHx(2:obj.Nx,:) = 2*Kx(1:obj.Nx-1,:).*Kx(2:obj.Nx,:)./(Kx(1:obj.Nx-1,:)+Kx(2:obj.Nx,:));
            KHy(:,2:obj.Ny) = 2*Ky(:,1:obj.Ny-1).*Ky(:,2:obj.Ny)./(Ky(:,1:obj.Ny-1)+Ky(:,2:obj.Ny));
            
            %Transmissibility
            obj.Tx(2:obj.Nx,:) = obj.Ax./obj.dx.*KHx(2:obj.Nx,:);
            obj.Ty(:,2:obj.Ny) = obj.Ay./obj.dy.*KHy(:,2:obj.Ny);
        end
        function AddCoordinates(obj)
            % Add I, J coordinates to Grid
            obj.I = ones(obj.N, 1);
            obj.J = ones(obj.N, 1);
            Jindexes = 1:1:obj.Ny;
            for i=1:obj.Ny
                a = obj.Nx*(i-1)+1;
                obj.I(a:a+obj.Nx-1) = 1:1:obj.Nx;
                obj.J(a:a+obj.Nx-1) = Jindexes(i)*ones(obj.Nx,1);
            end
        end
    end
end