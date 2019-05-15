% Gravity model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
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
        function obj = gravity_model(DiscretizationModel, n_phases, n_frac)
            obj.RhoInt = cell(n_phases, n_frac + 1);
            obj.initialize_rhoint(DiscretizationModel.ReservoirGrid, n_phases, 0);
            for f=1:n_frac
                obj.initialize_rhoint(DiscretizationModel.FracturesGrid.Grids(f), n_phases, f);
            end
        end
        function initialize_rhoint(obj, Grid, n_phases, f)
            nx = Grid.Nx;
            ny = Grid.Ny;
            nz = Grid.Nz;
            temp = struct('x', zeros(nx+1, ny, nz), 'y', zeros(nx, ny+1, nz), 'z', zeros(nx, ny, nz+1));
            for i=1:n_phases
                obj.RhoInt{i, f+1} = temp; 
            end
        end
        function ComputeInterfaceDensities(obj, Grid, Status, f)
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            Nz = Grid.Nz;
            [n_phases, ~] = size(obj.RhoInt);
            for i=1:n_phases
                rho1 = obj.g * reshape(Status.Properties(['rho_', num2str(i)]).Value, Nx, Ny, Nz);
                obj.RhoInt{i, f+1}.x(2:Nx,:,:) = (rho1(1:Nx-1,:,:) + rho1(2:Nx,:,:)) / 2;
                obj.RhoInt{i, f+1}.y(:,2:Ny,:) = (rho1(:,1:Ny-1,:) + rho1(:,2:Ny,:)) / 2;
                obj.RhoInt{i, f+1}.z(:,:,2:Nz) = (rho1(:,:,1:Nz-1) + rho1(:,:,2:Nz)) / 2;
            end
        end
    end
end

