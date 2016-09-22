% Gravity model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 20 September 2016
%Last modified: 22 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef gravity_model < handle
    properties (Constant)
        g = 9.806 % gravitational constant []
    end
    properties
        RhoInt1
        RhoInt2
    end
    methods
        function obj = gravity_model(nx, ny)
            obj.RhoInt1.x = zeros(nx+1, ny);
            obj.RhoInt1.y = zeros(nx, ny+1);
            obj.RhoInt2.x = zeros(nx+1, ny);
            obj.RhoInt2.y = zeros(nx, ny+1);
        end
        function ComputeInterfaceDensities(obj, Nx, Ny, rho)
            
            rho1 = obj.g * reshape(rho(:,2), Nx, Ny);
            obj.RhoInt1.x(2:obj.Nx,:) = (rho1(1:obj.Nx-1,:) + rho1(2:obj.Nx,:)) / 2;
            obj.RhoInt1.y(:,2:obj.Ny) = (rho1(:,1:obj.Ny-1) + rho1(:,2:obj.Ny)) / 2;
            
            rho2 = obj.g * reshape(rho(:,2), Nx, Ny);
            obj.RhoInt2.x(2:obj.Nx,:) = (rho2(1:obj.Nx-1,:) + rho2(2:obj.Nx,:)) / 2;
            obj.RhoInt2.y(:,2:obj.Ny) = (rho2(:,1:obj.Ny-1) + rho2(:,2:obj.Ny)) / 2;         
        end
    end
end

