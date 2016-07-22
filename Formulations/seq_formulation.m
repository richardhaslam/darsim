%  Sequential Formulation base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 22 July 2016
%Last modified: 22 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef seq_formulation < handle
    properties
        NumOfEquations
        Mob
        Mobt
        f
        U
    end
    methods
        function ComputeTotalMobility(obj, FluidModel, ProductionSystem)
            s = ProductionSystem.Reservoir.State.S;
            obj.Mob = FluidModel.ComputePhaseMobilities(s);
            obj.Mobt = sum(obj.Mob,2);
        end
        function UpdateFractionalFlow(obj, ProductionSystem, FluidModel)
            obj.f = obj.Mob(:,1) ./ sum(obj.Mob, 2);
        end
        function ComputeTransmissibilities(obj, DiscretizationModel, ProductionSystem)
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            dx = DiscretizationModel.ReservoirGrid.dx;
            dy = DiscretizationModel.ReservoirGrid.dy;
            Ax = DiscretizationModel.ReservoirGrid.Ax;
            Ay = DiscretizationModel.ReservoirGrid.Ay;
            K = ProductionSystem.Reservoir.K .* obj.Mobt;
            
            %Harmonic average of permeability.
            kx = reshape(K(:,1), Nx,Ny);
            ky = reshape(K(:,2), Nx, Ny);
            Kx = zeros(Nx, Ny+1);
            Ky = zeros(Nx, Ny+1);
            Kx(2:Nx,:) = 2*kx(1:Nx-1,:) .* kx(2:Nx,:) ./ (kx(1:Nx-1,:) + kx(2:Nx,:));
            Ky(:,2:Ny) = 2*ky(:,1:Ny-1) .* ky(:,2:Ny) ./ (ky(:,1:Ny-1) + ky(:,2:Ny));
            
            %Transmissibility
            obj.Tx = zeros(Nx+1, Ny);
            obj.Ty = zeros(Nx, Ny+1);
            obj.Tx(2:Nx,:) = Ax./dx .* Kx(2:Nx,:);
            obj.Ty(:,2:Ny) = Ay./dy .* Ky(:,2:Ny);
        end
        function BuildPressureResidual()
        end
        function BuildIncompressiblePressureMatrix(obj)
            
        end
        function ComputeFluxes(obj, ProductionSystem, DiscretizationModel)
            % 
            p = ProductionSystem.Reservoir.State.p;
            Nx = DiscretizationModel.ReservoirGrid.Nx;
            Ny = DiscretizationModel.ReservoirGrid.Ny;
            
            % Compute total fluxes ([m^3/s])
            P = reshape(p, Nx, Ny, 1);
            obj.U.x = zeros(Nx+1,Ny, 1);
            obj.U.y = zeros(Nx,Ny+1, 1);
            obj.U.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)) .* obj.Tx(2:Nx,:) - Ucap.x(2:Nx,:);
            obj.U.y(:,2:Ny) = (P(:,1:Ny-1)-P(:,2:Ny)) .* obj.Ty(:,2:Ny) - Ucap.y(:,2:Ny);
        end
        function BuildTransportResidual()
        end
        function UpdatePressure()
        end
        function UpdateSaturation()
        end
    end
end