%  Discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%  Discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Discretization_model < handle
    properties
        N
        ReservoirGrid
        FracturesGrid
        CrossConnections
    end
    methods
        function AddReservoirGrid(obj, reservoirgrid)
            obj.ReservoirGrid = reservoirgrid;
        end
        function AddFracturesGrid(obj, fracturesgrid)
            obj.FracturesGrid = fracturesgrid;
        end
        function AddCrossConnections(obj, crossconnections, Formulation)
            obj.CrossConnections = crossconnections; 
        end
        function Initialize(obj, ProductionSystem, Formulation)
            obj.ReservoirGrid.Initialize(ProductionSystem.Reservoir);    
            % Perforated cells
            obj.DefinePerforatedCells(ProductionSystem.Wells);
            
            % Total number of cells
            obj.N = obj.ReservoirGrid.N;
            
            % . Assign Depth
            obj.ReservoirGrid.ComputeDepth(Formulation.GravityModel.alpha, ProductionSystem.Reservoir.Thickness);
            obj.N = obj.ReservoirGrid.N;
            
            % initialize fractures
            if ProductionSystem.FracturesNetwork.Active
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    obj.FracturesGrid.Grids(f).Initialize(ProductionSystem.FracturesNetwork.Fractures(f));
                end
                % Total number of cells
                obj.N = obj.N + sum(obj.FracturesGrid.N);
                % Adding the harmonic permeabilities to CrossConnections
                obj.AddHarmonicPermeabilities(ProductionSystem.Reservoir, ProductionSystem.FracturesNetwork.Fractures);
            end            
        end
        function DefinePerforatedCells(obj, Wells)
            % Has to be improved for Diagonal wells (maybe using trajectories)
            % Injectors
            for i = 1:Wells.NofInj
                I = Wells.Inj(i).Coord(1,1):1:Wells.Inj(i).Coord(1,2);
                J = Wells.Inj(i).Coord(2,1):1:Wells.Inj(i).Coord(2,2);
                K = Wells.Inj(i).Coord(3,1):1:Wells.Inj(i).Coord(3,2);
                if (sum(I > obj.ReservoirGrid.Nx) == 1 || sum(J > obj.ReservoirGrid.Ny) == 1 || sum(K > obj.ReservoirGrid.Nz) == 1)
                    error(['ERROR: coordinates of injector num ', num2str(i),' fall outside of the domain']);
                else
                    Wells.Inj(i).Cells = I + (J-1)*obj.ReservoirGrid.Nx + (K-1)*obj.ReservoirGrid.Nx*obj.ReservoirGrid.Ny;
                    Wells.Inj(i).ResizeObjects(length(Wells.Inj(i).Cells));
                end
            end
            
            % Producers
            for i = 1:Wells.NofProd
                I = Wells.Prod(i).Coord(1,1):1:Wells.Prod(i).Coord(1,2);
                J = Wells.Prod(i).Coord(2,1):1:Wells.Prod(i).Coord(2,2);
                K = Wells.Prod(i).Coord(3,1):1:Wells.Prod(i).Coord(3,2);
                if (sum(I > obj.ReservoirGrid.Nx) == 1 || sum(J > obj.ReservoirGrid.Ny) == 1 || sum(K > obj.ReservoirGrid.Nz) == 1)
                    error(['ERROR: coordinates of producer num ', num2str(i),' fall outside of the domain']);
                else
                     Wells.Prod(i).Cells = I + (J-1)*obj.ReservoirGrid.Nx + (K-1)*obj.ReservoirGrid.Nx*obj.ReservoirGrid.Ny;
                     Wells.Prod(i).ResizeObjects(length(Wells.Prod(i).Cells));
                end
            end
        end
        function I = Index_Local_to_Global(obj, i, j, k, f, g)
            if (i<1),  error('i should at least be 1 but is not!');  end
            if (j<1),  error('j should at least be 1 but is not!');  end
            if (k<1),  error('k should at least be 1 but is not!');  end
            if (i>obj.ReservoirGrid.Nx),  i = obj.ReservoirGrid.Nx;  end
            if (j>obj.ReservoirGrid.Ny),  j = obj.ReservoirGrid.Ny;  end
            if (k>obj.ReservoirGrid.Nz),  k = obj.ReservoirGrid.Nz;  end
            if (f<0),  error('f (fracture index) cannot be negative!');  end
            if (f>obj.FracturesGrid.Nfrac),  error('f exceeds the number of fractures!');  end
            if (g<0),  error('g (fracture cell index) cannot be negative!');  end
            if (f~=0)&&(g>obj.FracturesGrid.Grids(f).N),  error('g exceeds the number of fracture cells!');  end
            
            if (f==0) && (g==0)
                I = (k-1)*(obj.ReservoirGrid.Nx*obj.ReservoirGrid.Ny) + (j-1)*(obj.ReservoirGrid.Nx) + i;
            elseif (f~=0) && (g~=0)
                I = (obj.ReservoirGrid.N) + sum(obj.FracturesGrid.N(1:f-1)) + g;
            else
                error('For global indexing in the reservoir both f & g must be zero!\nFor global indexing in the fractures both f & g must be non-zero!');
            end
        end
        function indexing = Index_Global_to_Local(obj, I)
            if (I<1),  error('Global indexing (I) cannot be negative!');  end
            if (I>obj.N),  error('Global indexing (I) cannot exceed total number of cells!');  end
            if I <= obj.ReservoirGrid.Nx*obj.ReservoirGrid.Ny*obj.ReservoirGrid.Nz
                indexing.i = mod( I , obj.ReservoirGrid.Nx );
                if ( indexing.i==0 ),   indexing.i = obj.ReservoirGrid.Nx;   end
                indexing.j = mod(  (I-indexing.i)/obj.ReservoirGrid.Nx  , obj.ReservoirGrid.Ny ) +1;
                if ( indexing.j==0 ),   indexing.j = obj.ReservoirGrid.Ny;   end
                indexing.k = mod( ((I-indexing.i)/obj.ReservoirGrid.Nx -indexing.j+1)/obj.ReservoirGrid.Ny , obj.ReservoirGrid.Nz ) +1;               
                if ( indexing.k==0 ),   indexing.k = obj.ReservoirGrid.Nz;   end
                indexing.f = 0;
                indexing.g = 0;
            else
                indexing.i = obj.ReservoirGrid.Nx;
                indexing.j = obj.ReservoirGrid.Ny;
                indexing.k = obj.ReservoirGrid.Nz;
                temp = I - obj.ReservoirGrid.N;
                temp = find( temp - cumsum(obj.FracturesGrid.N) <= 0);
                indexing.f = temp(1);
                indexing.g = I - obj.ReservoirGrid.N - sum( obj.FracturesGrid.N(1:indexing.f-1) );
                if indexing.g==0,  indexing.g = obj.FracturesGrid.Grids(indexing.f).N;  end
            end
            if Index_Local_to_Global(obj, indexing.i, indexing.j, indexing.k, indexing.f, indexing.g) ~= I
                error('i,j,k are not correspondent with I. Check the formula again!');
            end
        end
        function AddHarmonicPermeabilities(obj, Reservoir, Fractures)
            for If1_Local = 1 : length(obj.CrossConnections)
                If1_Global = obj.ReservoirGrid.N+If1_Local; % Global index of this fracture cell;
                Ind_frac1_Local = obj.Index_Global_to_Local(If1_Global);
                
                indices_m = obj.CrossConnections(If1_Local).Cells( obj.CrossConnections(If1_Local).Cells <= obj.ReservoirGrid.N );
                obj.CrossConnections(If1_Local).T_Geo(1:length(indices_m)) = obj.CrossConnections(If1_Local).ConnIndex(1:length(indices_m)) .* ...
                    ( (obj.ReservoirGrid.dx + obj.ReservoirGrid.dy + obj.ReservoirGrid.dz)/3 + Fractures(Ind_frac1_Local.f).Thickness ) ./...
                      ( ( (obj.ReservoirGrid.dx + obj.ReservoirGrid.dy + obj.ReservoirGrid.dz)/3 ./ Reservoir.K(indices_m,1) ) + ...
                        ( Fractures(Ind_frac1_Local.f).Thickness ./ Fractures(Ind_frac1_Local.f).K(Ind_frac1_Local.g,1) ) );
            
                indices_f = obj.CrossConnections(If1_Local).Cells( obj.CrossConnections(If1_Local).Cells > obj.ReservoirGrid.N );
                if ~isempty(indices_f)
                    for n = 1:length(indices_f)
                        If2_Global = indices_f(n); % Global indices of the other fractures' cells if any
                        If2_Local = If2_Global - obj.ReservoirGrid.N;
                        Ind_frac2_Local = obj.Index_Global_to_Local(If2_Global);
                        
                        obj.CrossConnections(If1_Local).T_Geo(length(indices_m)+n) = obj.CrossConnections(If1_Local).ConnIndex(length(indices_m)+n) * ...
                            ( (obj.FracturesGrid.Grids(Ind_frac1_Local.f).dx + obj.FracturesGrid.Grids(Ind_frac1_Local.f).dy)/2 + Fractures(Ind_frac2_Local.f).Thickness ) ./...
                              ( ( (obj.FracturesGrid.Grids(Ind_frac1_Local.f).dx + obj.FracturesGrid.Grids(Ind_frac1_Local.f).dy)/2 ./ Fractures(Ind_frac1_Local.f).K(Ind_frac1_Local.g,1) ) + ...
                                ( Fractures(Ind_frac2_Local.f).Thickness ./ Fractures(Ind_frac2_Local.f).K(Ind_frac2_Local.g,1) ) );
                    end
                end
            end
        end
        function AverageMassOnCoarseBlocks(obj, Status, FluidModel, Formulation)
            % virtual call
        end
    end
end