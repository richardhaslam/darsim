%  Grid mapper 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 5 September 2016
%Last modified: 5 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef grid_mapper < handle
    properties
    end
    methods
        function BuildFamily(obj, CoarseGrid, FineGrid, CF, level) 
        % Assign Children and GrandChildren 
          for c = 1:CoarseGrid.N
              %coordinates of fine cells contained in the coarse block
              J = ceil(c/CoarseGrid.Ny);
              I = c - CoarseGrid.Nx*(J - 1);
              Imin = (I - 1) * CF(2) + 1;
              Imax = Imin + CF(1)-1;
              Jmin = (J - 1) * CF(1) + 1;
              Jmax = Jmin + CF(2)-1;
              %indexes of the fine cells
              i = Imin:Imax;
              j = Jmin:Jmax;
              [p,q] = meshgrid(i, j);
              pairs = [p(:), q(:)];
              indexes = pairs(:,1) + (pairs(:,2)-1)*FineGrid.Nx;
              CoarseGrid.Children(c, :) = indexes;
              FineGrid.Fathers(indexes) = c;
              
              %coordinates of fine cells contained in the coarse block
              Imin = CoarseGrid.I(c) - floor((CF(1)^level - 1)/2);
              Imax = CoarseGrid.I(c) + ceil((CF(1)^level - 1)/2);
              Jmin = CoarseGrid.J(c) - floor((CF(2)^level - 1)/2);
              Jmax = CoarseGrid.J(c) + ceil((CF(2)^level - 1)/2);
              i = Imin:Imax;
              j = Jmin:Jmax;
              [p,q] = meshgrid(i, j);
              pairs = [p(:), q(:)];
              %indexes of the fine cells
              indexes = pairs(:,1) + (pairs(:,2)-1)*CoarseGrid.Nx*CF(1)^level;
              CoarseGrid.GrandChildren(c, :) = indexes;
          end          
        end
    end
end
    