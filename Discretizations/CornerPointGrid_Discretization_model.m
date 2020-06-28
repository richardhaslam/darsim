% CornerPointGrid Discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef CornerPointGrid_Discretization_model < FS_Discretization_model
    properties
        CornerPointGridData
    end
    methods
        function Initialize(obj, ProductionSystem, Formulation)
            obj.ReservoirGrid.Initialize(ProductionSystem.Reservoir);
            obj.ReservoirGrid.AddGridCoordinates();
            obj.ReservoirGrid.CorrectTransmissibilitiesForpEDFM();
            % Perforated cells
            obj.DefinePerforatedCells(ProductionSystem.Wells);
            
            % Total number of cells
            obj.N = obj.ReservoirGrid.N;
            
            % Assign Depth
            obj.ReservoirGrid.ComputeDepth(Formulation.GravityModel.alpha, ProductionSystem.Reservoir.Thickness);
            
            % Initialize fractures
            if ProductionSystem.FracturesNetwork.Active
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    obj.FracturesGrid.Grids(f).Initialize(ProductionSystem.FracturesNetwork.Fractures(f));
                    obj.FracturesGrid.Grids(f).CorrectTransmissibilitiesForpEDFM();
                end
                % Total number of cells
                obj.N = obj.N + sum(obj.FracturesGrid.N);
                % Adding the harmonic permeabilities to CrossConnections
                obj.AddHarmonicPermeabilities(ProductionSystem.Reservoir, ProductionSystem.FracturesNetwork.Fractures);
                % Adding the harmonic conductivities to CrossConnections
                if ~isempty(ProductionSystem.Reservoir.K_Cond_eff)
                    obj.AddHarmonicConductivities(ProductionSystem.Reservoir, ProductionSystem.FracturesNetwork.Fractures);
                end
            end
            
            % Merge the reservoir and all the fracture grids (if any) into one unified grid
            if ProductionSystem.FracturesNetwork.Active
                obj.AddUnifiedGrid();
            else
                obj.UnifiedGrid = obj.ReservoirGrid;
            end
        end
        function AddUnifiedGrid(obj)
            CPGData = obj.ReservoirGrid.CornerPointGridData;
            CPGData.Nx = [ obj.ReservoirGrid.Nx ; [obj.FracturesGrid.Grids.Nx]' ];
            CPGData.Ny = [ obj.ReservoirGrid.Ny ; [obj.FracturesGrid.Grids.Ny]' ];
            CPGData.Nz = [ obj.ReservoirGrid.Nz ; [obj.FracturesGrid.Grids.Nz]' ];
            
            Nf = obj.FracturesGrid.N;
            N_AllFaces_x = zeros( obj.FracturesGrid.Nfrac , 1);
            N_AllFaces_y = zeros( obj.FracturesGrid.Nfrac , 1);
            N_InternalFaces_x = zeros( obj.FracturesGrid.Nfrac , 1);
            N_InternalFaces_y = zeros( obj.FracturesGrid.Nfrac , 1);
            N_ExternalFaces_x = zeros( obj.FracturesGrid.Nfrac , 1);
            N_ExternalFaces_y = zeros( obj.FracturesGrid.Nfrac , 1);

            CPGData.N_ActiveCells = CPGData.N_ActiveCells + sum(Nf);

            for f = 1 : obj.FracturesGrid.Nfrac
                obj.FracturesGrid.Grids(f).Volume = obj.FracturesGrid.Grids(f).Volume .* ones(obj.FracturesGrid.Grids(f).N,1);
                Grid = obj.FracturesGrid.Grids(f);
                N_AllFaces_x(f) = length(Grid.Tx(:));
                N_AllFaces_y(f) = length(Grid.Ty(:));
                Internal_Tx = Grid.Tx( 2:Grid.Nx , 1:Grid.Ny );
                Internal_Ty = Grid.Ty( 1:Grid.Nx , 2:Grid.Ny );
                N_InternalFaces_x(f) = length(Internal_Tx(:));
                N_InternalFaces_y(f) = length(Internal_Ty(:));
                N_ExternalFaces_x(f) = N_AllFaces_x(f) - N_InternalFaces_x(f);
                N_ExternalFaces_y(f) = N_AllFaces_y(f) - N_InternalFaces_y(f);
            end
            
            CPGData.N_InternalFaces = CPGData.N_InternalFaces + sum(N_InternalFaces_x) + sum(N_InternalFaces_y);
            CPGData.N_ExternalFaces = CPGData.N_ExternalFaces + sum(N_ExternalFaces_x) + sum(N_ExternalFaces_y);

            Frac_Trans = cell(obj.FracturesGrid.Nfrac,1);
            for f = 1 : obj.FracturesGrid.Nfrac
                Grid = obj.FracturesGrid.Grids(f);
                % Merging the cell data
                for j = 1 : Grid.Ny
                    for i = 1 : Grid.Nx
                        CellInd = (j-1) * Grid.Nx + i;
                        CornerInd_SW = (j-1) * (Grid.Nx+1) + i;
                        CornerInd_SE = (j-1) * (Grid.Nx+1) + i+1;
                        CornerInd_NW = (j  ) * (Grid.Nx+1) + i;
                        CornerInd_NE = (j  ) * (Grid.Nx+1) + i+1;
                        CPGData.Cell.SW_Top_Corner(   obj.ReservoirGrid.N + sum(Nf(1:f-1))+CellInd , : ) = Grid.GridCoords( CornerInd_SW , : );
                        CPGData.Cell.SE_Top_Corner(   obj.ReservoirGrid.N + sum(Nf(1:f-1))+CellInd , : ) = Grid.GridCoords( CornerInd_SE , : );
                        CPGData.Cell.NW_Top_Corner(   obj.ReservoirGrid.N + sum(Nf(1:f-1))+CellInd , : ) = Grid.GridCoords( CornerInd_NW , : );
                        CPGData.Cell.NE_Top_Corner(   obj.ReservoirGrid.N + sum(Nf(1:f-1))+CellInd , : ) = Grid.GridCoords( CornerInd_NE , : );
                        CPGData.Cell.SW_Bot_Corner(   obj.ReservoirGrid.N + sum(Nf(1:f-1))+CellInd , : ) = Grid.GridCoords( CornerInd_SW , : );
                        CPGData.Cell.SE_Bot_Corner(   obj.ReservoirGrid.N + sum(Nf(1:f-1))+CellInd , : ) = Grid.GridCoords( CornerInd_SE , : );
                        CPGData.Cell.NW_Bot_Corner(   obj.ReservoirGrid.N + sum(Nf(1:f-1))+CellInd , : ) = Grid.GridCoords( CornerInd_NW , : );
                        CPGData.Cell.NE_Bot_Corner(   obj.ReservoirGrid.N + sum(Nf(1:f-1))+CellInd , : ) = Grid.GridCoords( CornerInd_NE , : );
                        CPGData.Cell.Centroid(        obj.ReservoirGrid.N + sum(Nf(1:f-1))+CellInd , : ) = Grid.Centroids( CellInd , : );
                        CPGData.Cell.Volume(          obj.ReservoirGrid.N + sum(Nf(1:f-1))+CellInd , : ) = Grid.Volume( CellInd );
                        CPGData.Cell.Index_Neighbors{ obj.ReservoirGrid.N + sum(Nf(1:f-1))+CellInd     } = Grid.Neighbours{ CellInd } + obj.ReservoirGrid.N + sum(Nf(1:f-1));
                        CPGData.Cell.N_Neighbors(     obj.ReservoirGrid.N + sum(Nf(1:f-1))+CellInd     ) = length( Grid.Neighbours{ CellInd } );
                    end
                end
                
                % Merging the internal Tx to the interfaces
                for j = 1 : Grid.Ny
                    for i = 1 : Grid.Nx-1
                        InterfaceInd = (j-1) * (Grid.Nx-1) + i;
                        CellNeighbor1Index = (j-1) * (Grid.Nx) + i;
                        CellNeighbor2Index = (j-1) * (Grid.Nx) + i+1;
                        
                        Frac_Trans{f} = vertcat( Frac_Trans{f} , Grid.Tx(i+1,j) );
                        
                        CPGData.Internal_Face.Area = vertcat( CPGData.Internal_Face.Area , Grid.dy * Grid.dz );
                        CPGData.Internal_Face.CellNeighbor1Index = vertcat( CPGData.Internal_Face.CellNeighbor1Index , CellNeighbor1Index + obj.ReservoirGrid.N + sum(Nf(1:f-1)) );
                        CPGData.Internal_Face.CellNeighbor2Index = vertcat( CPGData.Internal_Face.CellNeighbor2Index , CellNeighbor2Index + obj.ReservoirGrid.N + sum(Nf(1:f-1)) );
                    end
                end
                
                % Merging the internal Ty to the interfaces
                for i = 1 : Grid.Nx
                    for j = 1 : Grid.Ny-1
                        InterfaceInd = (i-1) * (Grid.Ny-1) + j;
                        CellNeighbor1Index = (j-1) * (Grid.Nx) + i;
                        CellNeighbor2Index = (j  ) * (Grid.Nx) + i;
                        
                        Frac_Trans{f} = vertcat( Frac_Trans{f} , Grid.Ty(i,j+1) );
                        
                        CPGData.Internal_Face.Area = vertcat( CPGData.Internal_Face.Area , Grid.dx * Grid.dz );
                        CPGData.Internal_Face.CellNeighbor1Index = vertcat( CPGData.Internal_Face.CellNeighbor1Index , CellNeighbor1Index + obj.ReservoirGrid.N + sum(Nf(1:f-1)) );
                        CPGData.Internal_Face.CellNeighbor2Index = vertcat( CPGData.Internal_Face.CellNeighbor2Index , CellNeighbor2Index + obj.ReservoirGrid.N + sum(Nf(1:f-1)) );
                    end
                end
            end
            
            
            % Creating the UnifiedGrid
            TempProperties.Grid.N = [obj.ReservoirGrid.Nx ; obj.ReservoirGrid.Ny ; obj.ReservoirGrid.Nz];
            TempProperties.CornerPointGridData = [];
            TempProperties.Grid.N_ActiveCells = obj.ReservoirGrid.N;
            TempProperties.CornerPointGridData.N_InternalFaces = length(obj.ReservoirGrid.Trans);
            obj.UnifiedGrid = corner_point_grid(TempProperties);     
            obj.UnifiedGrid.CornerPointGridData = CPGData;
            obj.UnifiedGrid.Nx = CPGData.Nx;
            obj.UnifiedGrid.Ny = CPGData.Ny;
            obj.UnifiedGrid.Nz = CPGData.Nz;
            obj.UnifiedGrid.N = sum(obj.UnifiedGrid.Nx .* obj.UnifiedGrid.Ny .* obj.UnifiedGrid.Nz);
            obj.UnifiedGrid.Volume = [ obj.ReservoirGrid.Volume ; vertcat(obj.FracturesGrid.Grids.Volume) ];
            obj.UnifiedGrid.Depth  = [ obj.ReservoirGrid.Depth  ; vertcat(obj.FracturesGrid.Grids.Depth ) ];
            obj.UnifiedGrid.Neighbours = obj.UnifiedGrid.CornerPointGridData.Cell.Index_Neighbors;
            obj.UnifiedGrid.Trans = [ obj.ReservoirGrid.Trans ; vertcat(Frac_Trans{:}) ];
            
            % Adding the non-neighboring CrosscConnections data to UnifiedGrid
            for If1Local = 1 : length(obj.CrossConnections)
                If1Global = obj.ReservoirGrid.N + If1Local;
                Cells = obj.CrossConnections(If1Local).Cells;
                T_Geo = obj.CrossConnections(If1Local).T_Geo;
                
                obj.UnifiedGrid.CornerPointGridData.N_InternalFaces = obj.UnifiedGrid.CornerPointGridData.N_InternalFaces + length(Cells);
                obj.UnifiedGrid.CornerPointGridData.Internal_Face.CellNeighbor1Index = vertcat( obj.UnifiedGrid.CornerPointGridData.Internal_Face.CellNeighbor1Index , If1Global*ones(length(Cells),1) );
                obj.UnifiedGrid.CornerPointGridData.Internal_Face.CellNeighbor2Index = vertcat( obj.UnifiedGrid.CornerPointGridData.Internal_Face.CellNeighbor2Index , Cells );
                obj.UnifiedGrid.Trans = vertcat( obj.UnifiedGrid.Trans , T_Geo );
            end
            
            obj.UnifiedGrid.ConstructConnectivityMatrix();
        end
        function AddHarmonicPermeabilities(obj, Reservoir, Fractures)

            Nm = obj.ReservoirGrid.N;
            for If1Local = 1 : length(obj.CrossConnections)
                Cells = obj.CrossConnections(If1Local).Cells;
                CI = obj.CrossConnections(If1Local).CI;
                If1Global = Nm+If1Local; % Global index of this fracture cell;
                Ind_frac1_Local = obj.Index_Global_to_Local(If1Global);
                f1 = Ind_frac1_Local.f;
                g1 = Ind_frac1_Local.g;
                
                indices_m = Cells( Cells <= Nm );
                obj.CrossConnections(If1Local).T_Geo(1:length(indices_m)) = ...
                    CI(1:length(indices_m)) .* ( obj.ReservoirGrid.Volume(indices_m).^(1/3) + Fractures(f1).Thickness ) ./...
                    ( ( obj.ReservoirGrid.Volume(indices_m).^(1/3) ./ Reservoir.K(indices_m,1) ) + ( Fractures(f1).Thickness ./ Fractures(f1).K(g1,1) ) );
   
                indices_f = Cells( Cells > Nm );
                if ~isempty(indices_f)
                    for n = 1:length(indices_f)
                        If2Global = indices_f(n);
                        If2Local = If2Global - Nm;
                        Ind_frac2_Local = obj.Index_Global_to_Local(indices_f(n));
                        f2 = Ind_frac2_Local.f;
                        g2 = Ind_frac2_Local.g;
                        obj.CrossConnections(If1Local).T_Geo(length(indices_m)+n) = ...
                            CI(length(indices_m)+n) * ( (obj.FracturesGrid.Grids(f1).dx + obj.FracturesGrid.Grids(f1).dy)/2 + Fractures(f2).Thickness ) ./...
                              ( ( (obj.FracturesGrid.Grids(f1).dx + obj.FracturesGrid.Grids(f1).dy)/2 ./ Fractures(f1).K(g1,1) ) + ( Fractures(f2).Thickness ./ Fractures(f2).K(g2,1) ) );
                    end
                end
            end
        end
%         function indexing = Index_Global_to_Local(obj, I)
%             if (I<1),  error('Global indexing (I) cannot be negative!');  end
%             if (I>obj.N),  error('Global indexing (I) cannot exceed total number of cells!');  end
%             if I <= obj.ReservoirGrid.N
%                 indexing.Im = I;
%                 indexing.f = 0;
%                 indexing.g = 0;
%             else
%                 indexing.Im = obj.ReservoirGrid.N;
%                 temp = I - obj.ReservoirGrid.N;
%                 temp = find( temp - cumsum(obj.FracturesGrid.N) <= 0);
%                 indexing.f = temp(1);
%                 indexing.g = I - obj.ReservoirGrid.N - sum( obj.FracturesGrid.N(1:indexing.f-1) );
%                 if indexing.g==0,  indexing.g = obj.FracturesGrid.Grids(indexing.f).N;  end
%             end
%             if obj.Index_Local_to_Global(indexing.Im, indexing.f, indexing.g) ~= I
%                 error('Im is not correspondent with I. Check the formula again!');
%             end
%         end
        function ObtainPerforatedCellsBasedOnIJKList(obj, Well, Well_Type, w)
            CellListIndex = [];
            for p = 1:size(Well.Coordinate.Value,1) - 1
                I = Well.Coordinate.Value(p,1) : 1 : Well.Coordinate.Value(p+1,1);
                J = Well.Coordinate.Value(p,2) : 1 : Well.Coordinate.Value(p+1,2);
                K = Well.Coordinate.Value(p,3) : 1 : Well.Coordinate.Value(p+1,3);
                if (sum(I > obj.ReservoirGrid.Nx) == 1 || sum(J > obj.ReservoirGrid.Ny) == 1 || sum(K > obj.ReservoirGrid.Nz) == 1)
                    error('DARSim Error: The coordinates of the %s #%d fall outside of the domain. Check the input file!\n', Well_Type, w);
                else
                CellListIndex = [ CellListIndex, I + (J-1)*obj.ReservoirGrid.Nx + (K-1)*obj.ReservoirGrid.Nx*obj.ReservoirGrid.Ny ];
                end
            end
            Well.Cells = zeros(length(CellListIndex),1);
            Well.Cells(:,1) = CellListIndex;
            Well.Cells = unique(Well.Cells);
            Well.ResizeObjects(length(Well.Cells));
        end
        function ObtainPerforatedCellsBasedOnXYZList(obj, Well, Well_Type, w)
            for p = 1:size(Well.Coordinate.Value,1) - 1
                PointA = Well.Coordinate.Value(p  ,:)';
                PointB = Well.Coordinate.Value(p+1,:)';
                LineSegment = lineSegment_DARSim(PointA,PointB);
                
                [ ~ , indList ] = min( vecnorm(LineSegment.PointM' - obj.CornerPointGridData.Cell.Centroid, 2,2) );
                Count = 1;
                while Count <= length(indList)
                    I = indList(Count);
                    Count = Count+1;
                    NW_Top = obj.CornerPointGridData.Cell.NW_Top_Corner(I,:)';
                    SW_Top = obj.CornerPointGridData.Cell.SW_Top_Corner(I,:)';
                    SE_Top = obj.CornerPointGridData.Cell.SE_Top_Corner(I,:)';
                    NE_Top = obj.CornerPointGridData.Cell.NE_Top_Corner(I,:)';
                    NW_Bot = obj.CornerPointGridData.Cell.NW_Bot_Corner(I,:)';
                    SW_Bot = obj.CornerPointGridData.Cell.SW_Bot_Corner(I,:)';
                    SE_Bot = obj.CornerPointGridData.Cell.SE_Bot_Corner(I,:)';
                    NE_Bot = obj.CornerPointGridData.Cell.NE_Bot_Corner(I,:)';
                    
                    Hexahedron = hexahedron_DARSim(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot);
                    Hexahedron.Centroid = obj.CornerPointGridData.Cell.Centroid(I,:)';
                    
                    Epsilon = 1e-10 * ( min(obj.ReservoirGrid.Volume) )^(1/3);
                    [Geostatus, IntersectPoints] = Hexahedron.Obtain_Polyhedron_LineSegment_Intersection(LineSegment, Epsilon);
                    
                    % Add the neighboring cells to the list for intersection check if it is the first try
                    % but no intersection occurs, so another one must be checked
                    if isempty(IntersectPoints) && Count == 2
                        indNeighbors = obj.CornerPointGridData.Cell.Index_Neighbors{I};
                        indList = union(indList, indNeighbors, 'stable');
                        continue;
                    end
                    
                    % Assigning the index of cell to the array
                    if isempty(IntersectPoints) && Count > 2
                        continue;
                    end
                    indNeighbors = obj.CornerPointGridData.Cell.Index_Neighbors{I};
                    indList = union(indList, indNeighbors, 'stable');
                    
                    if Geostatus.haveIntersect
                        Well.Cells = [ Well.Cells ; I ];
                    end
                end
            end
            Well.Cells = unique(Well.Cells);
            Well.ResizeObjects(length(Well.Cells));
        end
        function ObtainPerforatedCellsBasedOnCellIndexList(obj, Well, Well_Type, w)
            Well.Cells = zeros(length(Well.Coordinate.Value),1);
            Well.Cells(:,1) = Well.Coordinate.Value;
            Well.Cells = unique(Well.Cells);
            Well.ResizeObjects(length(Well.Cells));
        end
    end
end