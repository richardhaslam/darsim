%  MS discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Multiscale_Discretization_model < Discretization_model
    properties
        Nf
        Nc
        Vertex_On_Corner
        maxLevel
        CoarseningSwitch % This is a switch for Coarsening options (mostly valid for fractures).
        % "0" means coarsened upto the max coarsening level, afterwards stays at that resolution (higher computational demand)
        % "1" means coarsened upto the max coarsening level, afterwards no coare node (lower computational demand)
        Coarsening
        CoarseGrid
        OperatorsHandler
        FineGrid
        GridMapper
    end
    methods
        function obj = Multiscale_Discretization_model(maxlevel, coarsening)
            n = length(maxlevel);
            obj.CoarseGrid = coarse_grid.empty;
            obj.Coarsening = coarsening;
            obj.maxLevel   = maxlevel;
            obj.Nf = zeros(n, 1);
            obj.Nc = zeros(n, maxlevel(1));
            obj.GridMapper = grid_mapper();
        end
        function AddOperatorsHandler(obj, operatorshandler)
            obj.OperatorsHandler = operatorshandler;
        end
        function InitializeMapping(obj, ProductionSystem, FluidModel)
            disp([num2str(obj.maxLevel(1)), ' levels multiscale run']);
            % Construct Coarse Grids
            disp(newline);
            disp('Constructing coarse grids');
            obj.ConstructCoarseGrids(ProductionSystem.Wells.Inj, ProductionSystem.Wells.Prod);
            
            if ProductionSystem.FracturesNetwork.Active
                obj.Nf = [obj.ReservoirGrid.N; obj.FracturesGrid.N];
                obj.FineGrid = [obj.ReservoirGrid, obj.FracturesGrid.Grids];
            else
                obj.FineGrid = obj.ReservoirGrid;
                obj.Nf = obj.ReservoirGrid.N;
            end

            %% Pressure interpolators
            disp('Multiscale basis functions - start computation');
            start = tic;
            obj.OperatorsHandler.BuildMMsOperators(ProductionSystem, FluidModel, obj.FineGrid, obj.CrossConnections, ...
                obj.maxLevel, obj.CoarseGrid);
            disp('Multiscale basis functions - end')
            timer = toc(start);
            disp(['Multiscale basis functions computation took ', num2str(timer)])
            disp(newline);
        end
        function ConstructCoarseGrids(obj, Inj, Prod)
            obj.CheckIfNumberOfGridCellsMatchTheCoarseningSettings();
            
            % Maximum coarsening level is defined by resevoir. Fractures cannot have a higher coarsening level than reservoir:
            obj.maxLevel( obj.maxLevel > obj.maxLevel(1) ) = obj.maxLevel(1);
            
            %% 1. Reservoir
            % Construct all coarse grids for reservoir
            obj.ReservoirGrid.CoarseLevel = 0;
            obj.CoarseGrid(1,1) = coarse_grid();
            obj.CoarseGrid(1,1).CoarseLevel = 1;
            obj.CoarseGrid(1,1).hasCoarseNodes = 1;
            obj.CoarseGrid(1,1).CoarseFactor = obj.Coarsening(1,:,1);
            obj.CoarseGrid(1,1).Vertex_On_Corner = obj.Vertex_On_Corner;
            obj.CoarseGrid(1,1).BuildCoarseGrid(obj.ReservoirGrid, obj.Coarsening(1,:,1));
            obj.GridMapper.BuildFamily(obj.CoarseGrid(1,1), obj.ReservoirGrid, obj.Coarsening(1,:,1), 1);
            obj.CoarseGrid(1,1).AddWells(Inj, Prod);
            obj.Nc(1, 1) = obj.CoarseGrid(1,1).N;
            for L=2:obj.maxLevel(1)
                obj.CoarseGrid(1,L) = coarse_grid();
                obj.CoarseGrid(1,L).CoarseLevel = L;
                obj.CoarseGrid(1,L).hasCoarseNodes = 1;
                obj.CoarseGrid(1,L).CoarseFactor = obj.Coarsening(1,:,L);
                obj.CoarseGrid(1,L).Vertex_On_Corner = obj.Vertex_On_Corner;
                obj.CoarseGrid(1,L).BuildCoarseGrid([obj.ReservoirGrid, obj.CoarseGrid(1,1:L-1)],obj.Coarsening(1,:,1:L));
                obj.GridMapper.BuildFamily(obj.CoarseGrid(1,L), [obj.ReservoirGrid, obj.CoarseGrid(1,1:L-1)], obj.Coarsening(1,:,L)./obj.Coarsening(1,:,L-1), L);
                obj.CoarseGrid(1,L).AddWells(Inj, Prod);
                obj.Nc(1, L) = obj.CoarseGrid(1,L).N;
            end
            
            % Fathers and Verteces for reservoir
            obj.GridMapper.AssignFathersandVerteces([obj.ReservoirGrid, obj.CoarseGrid(1,:)], obj.maxLevel(1))
            
            %% 2. Fractures
            % Construct all coarse grids for fractures
            for f = 1 : size(obj.Coarsening,1) - 1
                % minMaxLevel = min( obj.maxLevel(1) , obj.maxLevel(1+f) );
                obj.FracturesGrid.Grids(f).CoarseLevel = 0;
                obj.CoarseGrid(1+f,1) = coarse_grid();
                obj.CoarseGrid(1+f,1).CoarseLevel = 1;
                obj.CoarseGrid(1+f,1).hasCoarseNodes = 1;
                obj.CoarseGrid(1+f,1).CoarseFactor = obj.Coarsening(1+f,:,1);
                obj.CoarseGrid(1+f,1).Vertex_On_Corner = obj.Vertex_On_Corner;
                obj.CoarseGrid(1+f,1).BuildCoarseGrid(obj.FracturesGrid.Grids(f),obj.Coarsening(1+f,:,1));
                obj.GridMapper.BuildFamily(obj.CoarseGrid(1+f,1), obj.FracturesGrid.Grids(f), obj.Coarsening(1+f,:,1), 1);
                obj.Nc(f+1,1) = obj.CoarseGrid(1+f,1).N;
                for L = 2 : obj.maxLevel(1)
                    obj.CoarseGrid(1+f,L) = coarse_grid();
                    obj.CoarseGrid(1+f,L).CoarseLevel = L;
                    obj.CoarseGrid(1+f,L).hasCoarseNodes = 1;
                    obj.CoarseGrid(1+f,L).CoarseFactor = obj.Coarsening(1+f,:,L);
                    obj.CoarseGrid(1+f,L).Vertex_On_Corner = obj.Vertex_On_Corner;
                    obj.CoarseGrid(1+f,L).BuildCoarseGrid([obj.FracturesGrid.Grids(f), obj.CoarseGrid(1+f,1:L-1)],obj.Coarsening(1+f,:,1:L));
                    obj.GridMapper.BuildFamily(obj.CoarseGrid(1+f,L), [obj.FracturesGrid.Grids(f), obj.CoarseGrid(1+f,1:L-1)], obj.Coarsening(1+f,:,L)./obj.Coarsening(1+f,:,L-1), L);
                    obj.Nc(f+1,L) = obj.CoarseGrid(1+f,L).N;
                end
                obj.GridMapper.AssignFathersandVerteces([obj.FracturesGrid.Grids(f), obj.CoarseGrid(1+f,:)], obj.maxLevel(1));
            end
            
            %% 3. Adding the non-neighboring connections to the list of neighbors for each coarse nodes (only for fractured cases)
            if size(obj.CoarseGrid,1) > 1
                obj.AddNonNeighboringConnectionsToNeighborsListAtCoarseScale();
            end
            
            %% 4. Flagging the coarse grids of fractures whose coarse nodes will not be used in ADM
            for f = 1 : size(obj.Coarsening,1) - 1
                if obj.CoarseningSwitch(1+f) == 1
                    % This fracture is coarsened upto the last coarsening level,
                    % afterwards it will disappear (lowest computational demand).
                    for L = obj.maxLevel(1+f)+1 : obj.maxLevel(1)
                        obj.CoarseGrid(1+f,L).hasCoarseNodes = 0;
                        obj.CoarseGrid(1+f,L).Verteces(:) = 0;
                    end
                    obj.FracturesGrid.Grids(f).Verteces( : , obj.maxLevel(1+f)+1:end ) = 0;
                    for L = 1 : obj.maxLevel(1+f)
                        obj.CoarseGrid(1+f,L).Verteces( : , obj.maxLevel(1+f)-L+1:end ) = 0;
                    end
                end
            end

            %% 4. Printing out the coarsening data
			fprintf('Coarsening ratio in reservoir: %d x %d x %d\n' , obj.Coarsening(1,1,1), obj.Coarsening(1,2,1), obj.Coarsening(1,3,1) );
            for L = 1 : obj.maxLevel(1)
                fprintf('Number of reservoir coarse nodes at level %d: %d\n' , L, obj.Nc(1,L).*obj.CoarseGrid(1,L).hasCoarseNodes );
                if (size(obj.Coarsening,1) - 1) > 0
                    fprintf('Number of fractures coarse nodes at level %d: %d\n' , L, sum(obj.Nc(2:end,L).*[obj.CoarseGrid(2:end,L).hasCoarseNodes]') );
                end
            end
        end
        function AddWellsToInitialPressure(obj, ProductionSystem, FluidModel)
            % Improving the first pressure guess for multilevel/multiscale method
            deltaP_w = obj.OperatorsHandler.ProlongationBuilders(1).StaticMultilevelPressureGuess(ProductionSystem, FluidModel, obj.FineGrid, obj.CoarseGrid(:, end), obj.CrossConnections);

            %% 1. Updating Reservoir Pressure with Well Correction
            Pm = ProductionSystem.Reservoir.State.Properties(['P_', num2str(FluidModel.NofPhases)]);
            Pm.update(deltaP_w(1:obj.ReservoirGrid.N));
            
            %% 2. Update fractures pressure and densities
            if ProductionSystem.FracturesNetwork.Active
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    % Update Pressure
                    index1 = obj.ReservoirGrid.N + sum(obj.FracturesGrid.N(1:f-1)) + 1;
                    index2 = index1 - 1 + obj.FracturesGrid.N(f);
                    Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_', num2str(FluidModel.NofPhases)]);
                    Pf.update(deltaP_w(index1:index2));
                end
            end
            
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
        end
        function SelectADMGrid(obj, ProductionSystem)
            % virtual call
        end
        function CheckIfNumberOfGridCellsMatchTheCoarseningSettings(obj)
            % Check for reservoir
            if mod( obj.ReservoirGrid.Nx , obj.Coarsening(1,1,end) ) == 0 && ...
               mod( obj.ReservoirGrid.Ny , obj.Coarsening(1,2,end) ) == 0 && ...
               mod( obj.ReservoirGrid.Nz , obj.Coarsening(1,3,end) ) == 0
                obj.Vertex_On_Corner = 0;
                fprintf('The verteces are not on the corners.\n');
            elseif mod( max((obj.ReservoirGrid.Nx-1),1) , obj.Coarsening(1,1,end) ) == 0 && ...
                   mod( max((obj.ReservoirGrid.Ny-1),1) , obj.Coarsening(1,2,end) ) == 0 && ...
                   mod( max((obj.ReservoirGrid.Nz-1),1) , obj.Coarsening(1,3,end) ) == 0
                obj.Vertex_On_Corner = 1;
                fprintf('The verteces are on the corners.\n');
            else
                error('DARSim2 Error: The number of reservoir grid cells do not comply with coarsening settings. Please check the input files.');
            end
            % Check for fractures
            for f = 1 : size(obj.Coarsening,1) - 1
                if mod( obj.FracturesGrid.Grids(f).Nx , obj.Coarsening(1+f,1,end) ) == 0 && ...
                   mod( obj.FracturesGrid.Grids(f).Ny , obj.Coarsening(1+f,2,end) ) == 0 && ...
                   mod( obj.FracturesGrid.Grids(f).Nz , obj.Coarsening(1+f,3,end) ) == 0
                   % It is OK. Do nothing.
                elseif mod( max((obj.FracturesGrid.Grids(f).Nx-1),1) , obj.Coarsening(1+f,1,end) ) == 0 && ...
                       mod( max((obj.FracturesGrid.Grids(f).Ny-1),1) , obj.Coarsening(1+f,2,end) ) == 0 && ...
                       mod( max((obj.FracturesGrid.Grids(f).Nz-1),1) , obj.Coarsening(1+f,3,end) ) == 0
                       % It is OK. Do nothing.
                else
                    error('DARSim2 Error: The number of grid cells in fracture #%d  do not comply with coarsening settings. Please check the input files.',f);
                end
            end
        end
        function AddNonNeighboringConnectionsToNeighborsListAtCoarseScale(obj)
            %% Initializing the "NonNeighbor" properties for all CoarseGrids
            for L = 1 : size(obj.CoarseGrid,2)
                for m = 1 : size(obj.CoarseGrid,1)
                    obj.CoarseGrid(m,L).NonNeighbours = cell(obj.CoarseGrid(m,L).N,1);
                end
            end

            %% 3.1 Coarsening level 1
            Nfx = [obj.ReservoirGrid.Nx; [obj.FracturesGrid.Grids.Nx]' ];
            Nfy = [obj.ReservoirGrid.Ny; [obj.FracturesGrid.Grids.Ny]' ];
            Nfz = [obj.ReservoirGrid.Nz; [obj.FracturesGrid.Grids.Nz]' ];
            Nf  = [obj.ReservoirGrid.N; obj.FracturesGrid.N];
            Ncx = [obj.CoarseGrid(:,1).Nx]';
            Ncy = [obj.CoarseGrid(:,1).Ny]';
            Ncz = [obj.CoarseGrid(:,1).Nz]';
            Nc  = [obj.CoarseGrid(:,1).N ]';
            
            % 3.1.1 Reservoir
            for c = 1 : obj.CoarseGrid(1,1).N
                Children = obj.CoarseGrid(1,1).Children{c,1};
                FineScaleNonNeighbors = unique([obj.ReservoirGrid.NonNeighbours{Children}]);
                for n = 1 : length( FineScaleNonNeighbors )
                    FineScaleLocalInd = obj.Index_Global_to_Local_CoarseScale( Nfx, Nfy, Nfz, FineScaleNonNeighbors(n) );
                    f = FineScaleLocalInd.f;  g = FineScaleLocalInd.g;
                    G = obj.FracturesGrid.Grids(f).Fathers(g,1);
                    CoarseScaleNonNeighbor = obj.Index_Local_to_Global_CoarseScale( Ncx, Ncy, Ncz, Ncx(1), Ncy(1), Ncz(1), f, G);
                    obj.CoarseGrid(1,1).NonNeighbours{c} = [obj.CoarseGrid(1,1).NonNeighbours{c}, CoarseScaleNonNeighbor];
                end
                obj.CoarseGrid(1,1).NonNeighbours{c} = unique( obj.CoarseGrid(1,1).NonNeighbours{c} );
            end
            
            % 3.1.2 Fractures
            for f = 1 : obj.FracturesGrid.Nfrac
                for c = 1 : obj.CoarseGrid(1+f,1).N
                    Children = obj.CoarseGrid(1+f,1).Children{c,1};
                    FineScaleNonNeighbors = unique([obj.FracturesGrid.Grids(f).NonNeighbours{Children}]);
                    % 3.1.2.1 non-neighboring cells from reservoir for this fracture
                    FineScaleNonNeighborsFromReservoir = FineScaleNonNeighbors( FineScaleNonNeighbors <= Nf(1) );
                    CoarseScaleNonNeighbor = obj.ReservoirGrid.Fathers(FineScaleNonNeighborsFromReservoir,1);
                    obj.CoarseGrid(1+f,1).NonNeighbours{c} = [ obj.CoarseGrid(1+f,1).NonNeighbours{c}, CoarseScaleNonNeighbor'];

                    % 3.1.2.2 non-neighboring cells from other fractures for this fracture
                    FineScaleNonNeighborsFromFractures = FineScaleNonNeighbors( FineScaleNonNeighbors > Nf(1) );
                    for n = 1 : length( FineScaleNonNeighborsFromFractures )
                        FineScaleLocalInd = obj.Index_Global_to_Local_CoarseScale( Nfx, Nfy, Nfz, FineScaleNonNeighborsFromFractures(n) );
                        f2 = FineScaleLocalInd.f;  g2 = FineScaleLocalInd.g;
                        G2 = obj.FracturesGrid.Grids(f2).Fathers(g2,1);
                        CoarseScaleNonNeighbor = obj.Index_Local_to_Global_CoarseScale( Ncx, Ncy, Ncz, Ncx(1), Ncy(1), Ncz(1), f2, G2);
                        obj.CoarseGrid(1+f,1).NonNeighbours{c} = [obj.CoarseGrid(1+f,1).NonNeighbours{c}, CoarseScaleNonNeighbor];
                    end
                    
                    obj.CoarseGrid(1+f,1).NonNeighbours{c} = unique( obj.CoarseGrid(1+f,1).NonNeighbours{c} );
                end
            end
            
            %% 3.2 Coarsening levels above 1
            for L = 2 : obj.maxLevel(1)
                Nfx = [obj.CoarseGrid(:,L-1).Nx]';
                Nfy = [obj.CoarseGrid(:,L-1).Ny]';
                Nfz = [obj.CoarseGrid(:,L-1).Nz]';
                Nf  = [obj.CoarseGrid(:,L-1).N ]';
                Ncx = [obj.CoarseGrid(:,L).Nx]';
                Ncy = [obj.CoarseGrid(:,L).Ny]';
                Ncz = [obj.CoarseGrid(:,L).Nz]';
                Nc  = [obj.CoarseGrid(:,L).N ]';
                
                % 3.2.1 reservoir
                for c = 1 : obj.CoarseGrid(1,L).N
                    Children = obj.CoarseGrid(1,L).Children{c,1};
                    FineScaleNonNeighbors = unique([obj.CoarseGrid(1,L-1).NonNeighbours{Children}]);
                    for n = 1 : length( FineScaleNonNeighbors )
                        FineScaleLocalInd = obj.Index_Global_to_Local_CoarseScale( Nfx, Nfy, Nfz, FineScaleNonNeighbors(n) );
                        f = FineScaleLocalInd.f;  g = FineScaleLocalInd.g;
                        G = obj.CoarseGrid(1+f,L-1).Fathers(g,1);
                        CoarseScaleNonNeighbor = obj.Index_Local_to_Global_CoarseScale( Ncx, Ncy, Ncz, Ncx(1), Ncy(1), Ncz(1), f, G);
                        obj.CoarseGrid(1,L).NonNeighbours{c} = [obj.CoarseGrid(1,L).NonNeighbours{c}, CoarseScaleNonNeighbor];
                    end
                    obj.CoarseGrid(1,L).NonNeighbours{c} = unique(obj.CoarseGrid(1,L).NonNeighbours{c});
                end
                
                % 3.2.2 Fractures
                for f = 1 : obj.FracturesGrid.Nfrac
                    for c = 1 : obj.CoarseGrid(1+f,L).N
                        Children = obj.CoarseGrid(1+f,L).Children{c,1};
                        FineScaleNonNeighbors = unique([obj.CoarseGrid(1+f,L-1).NonNeighbours{Children}]);
                        % 3.2.2.1 non-neighboring cells from reservoir for this fracture
                        FineScaleNonNeighborsFromReservoir = FineScaleNonNeighbors( FineScaleNonNeighbors <= Nf(1) );
                        CoarseScaleNonNeighbor = obj.CoarseGrid(1,L-1).Fathers(FineScaleNonNeighborsFromReservoir,1);
                        obj.CoarseGrid(1+f,L).NonNeighbours{c} = [ obj.CoarseGrid(1+f,L).NonNeighbours{c}, CoarseScaleNonNeighbor'];
                        
                        % 3.2.2.2 non-neighboring cells from other fractures for this fracture
                        FineScaleNonNeighborsFromFractures = FineScaleNonNeighbors( FineScaleNonNeighbors > Nf(1) );
                        for n = 1 : length( FineScaleNonNeighborsFromFractures )
                            FineScaleLocalInd = obj.Index_Global_to_Local_CoarseScale( Nfx, Nfy, Nfz, FineScaleNonNeighborsFromFractures(n) );
                            f2 = FineScaleLocalInd.f;  g2 = FineScaleLocalInd.g;
                            G2 = obj.CoarseGrid(1+f2,L-1).Fathers(g2,1);
                            CoarseScaleNonNeighbor = obj.Index_Local_to_Global_CoarseScale( Ncx, Ncy, Ncz, Ncx(1), Ncy(1), Ncz(1), f2, G2);
                            obj.CoarseGrid(1+f,L).NonNeighbours{c} = [obj.CoarseGrid(1+f,L).NonNeighbours{c}, CoarseScaleNonNeighbor];
                        end
                        
                        obj.CoarseGrid(1+f,L).NonNeighbours{c} = unique( obj.CoarseGrid(1+f,L).NonNeighbours{c} );
                    end
                end
            end
        end
        function I = Index_Local_to_Global_CoarseScale(obj, Nx, Ny, Nz, i, j, k, f, g)
            % In Nx,Ny,Nz, the first elements refers to reservoir and the rest are for fractures (1+f)
            N = Nx.*Ny.*Nz;
            if (i<1),  error('i should at least be 1 but is not!');  end
            if (j<1),  error('j should at least be 1 but is not!');  end
            if (k<1),  error('k should at least be 1 but is not!');  end
            if (i>Nx(1)),  i = Nx(1);  end
            if (j>Ny(1)),  j = Ny(1);  end
            if (k>Nz(1)),  k = Nz(1);  end
            if (f<0),  error('f (fracture index) cannot be negative!');  end
            if (f>length(N)-1),  error('f exceeds the number of fractures!');  end
            if (g<0),  error('g (fracture cell index) cannot be negative!');  end
            if (f~=0) && ( g > N(1+f) ),  error('g exceeds the number of fracture cells!');  end
            
            if (f==0) && (g==0)
                I = (k-1)*(Nx(1)*Ny(1)) + (j-1)*(Nx(1)) + i;
            elseif (f~=0) && (g~=0)
                I = sum(N(1:1+f-1)) + g;
            else
                error('For global indexing in the reservoir both f & g must be zero!\nFor global indexing in the fractures both f & g must be non-zero!');
            end
        end
        function indexing = Index_Global_to_Local_CoarseScale(obj, Nx, Ny, Nz, I)
            % In Nx,Ny,Nz, the first elements refers to reservoir and the rest are for fractures (1+f)
            N = Nx.*Ny.*Nz;
            if (I<1),  error('Global indexing (I) cannot be negative!');  end
            if (I>sum(N)),  error('Global indexing (I) cannot exceed total number of cells!');  end
            if I <= N(1)
                indexing.i = mod( I , Nx(1) );
                if ( indexing.i==0 ),   indexing.i = Nx(1);   end
                indexing.j = mod(  (I-indexing.i)/Nx(1) , Ny(1) ) +1;
                if ( indexing.j==0 ),   indexing.j = Ny(1);   end
                indexing.k = mod( ((I-indexing.i)/Nx(1) -indexing.j+1)/Ny(1) , Nz(1) ) +1;               
                if ( indexing.k==0 ),   indexing.k = Nz(1);   end
                indexing.f = 0;
                indexing.g = 0;
            else
                indexing.i = Nx(1);
                indexing.j = Ny(1);
                indexing.k = Nz(1);
                temp = I - N(1);
                temp = find( temp - cumsum(N(2:end)) <= 0);
                indexing.f = temp(1);
                indexing.g = I - N(1) - sum( N(2:indexing.f+1-1) );
                if indexing.g==0,  indexing.g = N(1+indexing.f);  end
            end
            if obj.Index_Local_to_Global_CoarseScale(Nx, Ny, Nz, indexing.i, indexing.j, indexing.k, indexing.f, indexing.g) ~= I
                error('i,j,k are not correspondent with I. Check the formula again!');
            end
        end
    end
end