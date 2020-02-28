% Gravity model class for CornerPointGrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: 
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef gravity_model_CornerPointGrid < gravity_model
    properties
    end
    methods
        function initialize_rhoint(obj, Grid, n_phases, f)
            for i=1:n_phases
                obj.RhoInt{i, f+1} = zeros(size(Grid.Trans)); 
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

