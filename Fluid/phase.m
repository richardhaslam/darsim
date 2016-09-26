% Phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 26 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef phase < handle
    properties
        rho0 % Reference density
        mu % Reference Viscosity
        cf % Compressibility
        sr % Irriducible saturation
    end
    properties (Constant)
        Pref = 1e5; % Atmospheric pressure
    end
    methods
        function [rho, mu] = UpdatePhaseProperties(obj, Status)
            rho = ComputeDensity(obj, Status);
            mu = obj.mu;
        end
        function rho = ComputeDensity(obj, p)
            rho = obj.rho0 .* exp(obj.cf.*(p - obj.Pref));
        end
        function [A, U] = UpWindAndRockFluxes(obj, Grid, P, RhoInt)
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            N = Grid.N;
            
            % Gravitational velocities
            depth = reshape(Grid.Depth, Nx, Ny);
            
            %
            Ugx = zeros(Nx+1,Ny,1);
            Ugy = zeros(Nx,Ny+1,1);
            Ugx(2:Nx,:) = (depth(1:Nx-1,:) - depth(2:Nx,:)) .* RhoInt.x(2:Nx,:); 
            Ugy(:,2:Ny) = (depth(:,1:Ny-1) - depth(:,2:Ny)) .* RhoInt.y(:,2:Ny);      
            
            P = reshape(P, Nx, Ny);
            %Compute 'rock' fluxes ([m^3/s])
            U.x = zeros(Nx+1,Ny,1);
            U.y = zeros(Nx,Ny+1,1);
            U.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)) .* Grid.Tx(2:Nx,:) + Grid.Tx(2:Nx,:) .* Ugx(2:Nx,:);
            U.y(:,2:Ny) = (P(:,1:Ny-1)-P(:,2:Ny)) .* Grid.Ty(:,2:Ny) + Grid.Ty(:,2:Ny) .* Ugy(:,2:Ny);
            
            %Use velocity to build upwind operator
            R = reshape((U.x(2:Nx+1,:) >= 0), N, 1);
            L = reshape((U.x(1:Nx,:) < 0), N, 1);
            T = reshape((U.y(:,2:Ny+1) >= 0), N, 1);
            B = reshape((U.y(:,1:Ny) < 0), N, 1);
            DiagVecs = [R, L];
            DiagIndx = [0,1];
            A.x = spdiags(DiagVecs,DiagIndx,N,N);
            DiagVecs = [T, B];
            DiagIndx = [0,Nx];
            A.y = spdiags(DiagVecs,DiagIndx,N,N);
        end
        function drho = DrhoDp(obj, p)
            drho = obj.cf .* obj.rho0 .*exp (obj.cf.*(p - obj.Pref));
        end
    end
end