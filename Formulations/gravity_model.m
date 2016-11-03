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
        function obj = gravity_model(nx, ny)
            obj.RhoInt = struct('x', zeros(nx+1, ny), 'y', zeros(nx, ny+1));
            obj.RhoInt(2).x = zeros(nx+1, ny);
            obj.RhoInt(2).y = zeros(nx, ny+1);
        end
        function ComputeInterfaceDensities(obj, Nx, Ny, rho)
            
            rho1 = obj.g * reshape(rho(:,1), Nx, Ny);
            obj.RhoInt(1).x(2:Nx,:) = (rho1(1:Nx-1,:) + rho1(2:Nx,:)) / 2;
            obj.RhoInt(1).y(:,2:Ny) = (rho1(:,1:Ny-1) + rho1(:,2:Ny)) / 2;
            
            rho2 = obj.g * reshape(rho(:,2), Nx, Ny);
            obj.RhoInt(2).x(2:Nx,:) = (rho2(1:Nx-1,:) + rho2(2:Nx,:)) / 2;
            obj.RhoInt(2).y(:,2:Ny) = (rho2(:,1:Ny-1) + rho2(:,2:Ny)) / 2;         
        end
    end
end

