% Gravity model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 20 September 2016
%Last modified: 22 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef gravity_model < handle
    properties
        g  % gravitational constant []
        alpha = 0; % in radiants
    end
    properties
        RhoInt
    end
    methods
        function obj = gravity_model(nx, ny, nz, n_phases)
            obj.RhoInt = struct('x', zeros(nx+1, ny, nz), 'y', zeros(nx, ny+1, nz), 'z', zeros(nx, ny, nz+1));
            for i=1:n_phases
                obj.RhoInt(i).x = zeros(nx+1, ny, nz);
                obj.RhoInt(i).y = zeros(nx, ny+1, nz);
                obj.RhoInt(i).z = zeros(nx, ny, nz+1);
            end
        end
        function ComputeInterfaceDensities(obj, Nx, Ny, Nz, Status)            
            for i=1:length(obj.RhoInt)
                rho1 = obj.g * reshape(Status.Properties(['rho_', num2str(i)]).Value, Nx, Ny, Nz);
                obj.RhoInt(i).x(2:Nx,:,:) = (rho1(1:Nx-1,:,:) + rho1(2:Nx,:,:)) / 2;
                obj.RhoInt(i).y(:,2:Ny,:) = (rho1(:,1:Ny-1,:) + rho1(:,2:Ny,:)) / 2;
                obj.RhoInt(i).z(:,:,2:Nz) = (rho1(:,:,1:Nz-1) + rho1(:,:,2:Nz)) / 2;
            end
        end
    end
end

