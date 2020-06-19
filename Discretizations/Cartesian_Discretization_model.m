% Cartesian Discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Cartesian_Discretization_model < FS_Discretization_model
    properties
    end
    methods
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
            Nx = obj.ReservoirGrid.Nx;
            Ny = obj.ReservoirGrid.Ny;
            Nz = obj.ReservoirGrid.Nz;
            dx = obj.ReservoirGrid.dx;
            dy = obj.ReservoirGrid.dy;
            dz = obj.ReservoirGrid.dz;
            Centroids = obj.ReservoirGrid.Centroids;
            
            for p = 1:size(Well.Coordinate.Value,1) - 1
                PointA = Well.Coordinate.Value(p  ,:)';
                PointB = Well.Coordinate.Value(p+1,:)';
                LineSegment = lineSegment_DARSim(PointA,PointB);
                [ ~ , indList ] = min( vecnorm(LineSegment.PointM' - Centroids, 2,2) );
                Count = 1;
                while Count <= length(indList)
                    I = indList(Count);
                    Count = Count+1;
                    NW_Top = [ Centroids(I,1) - dx/2 , Centroids(I,2) + dy/2 , Centroids(I,3) + dz/2 ]';
                    SW_Top = [ Centroids(I,1) - dx/2 , Centroids(I,2) - dy/2 , Centroids(I,3) + dz/2 ]';
                    SE_Top = [ Centroids(I,1) + dx/2 , Centroids(I,2) - dy/2 , Centroids(I,3) + dz/2 ]';
                    NE_Top = [ Centroids(I,1) + dx/2 , Centroids(I,2) + dy/2 , Centroids(I,3) + dz/2 ]';
                    NW_Bot = [ Centroids(I,1) - dx/2 , Centroids(I,2) + dy/2 , Centroids(I,3) - dz/2 ]';
                    SW_Bot = [ Centroids(I,1) - dx/2 , Centroids(I,2) - dy/2 , Centroids(I,3) - dz/2 ]';
                    SE_Bot = [ Centroids(I,1) + dx/2 , Centroids(I,2) - dy/2 , Centroids(I,3) - dz/2 ]';
                    NE_Bot = [ Centroids(I,1) + dx/2 , Centroids(I,2) + dy/2 , Centroids(I,3) - dz/2 ]';
                    
                    Hexahedron = hexahedron_DARSim(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot);
                    Hexahedron.Centroid = Centroids(I,:)';
                    
                    Epsilon = 1e-10 * ( min(obj.ReservoirGrid.Volume) )^(1/3);
                    [Geostatus, IntersectPoints] = Hexahedron.Obtain_Polyhedron_LineSegment_Intersection(LineSegment,Epsilon);
                    
                    % Add the neighboring cells to the list for intersection check if it is the first try
                    % but no intersection occurs, so another one must be checked
                    if isempty(IntersectPoints) && Count == 2
                        indNeighbors = obj.ReservoirGrid.Neighbours(I).indexes;
                        indList = union(indList, indNeighbors, 'stable');
                        continue;
                    end
                    
                    % Assigning the index of cell to the array
                    if isempty(IntersectPoints) && Count > 2
                        continue;
                    end
                    indNeighbors = obj.ReservoirGrid.Neighbours(I).indexes;
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
            Dxm = obj.ReservoirGrid.dx;
            Dym = obj.ReservoirGrid.dy;
            Dzm = obj.ReservoirGrid.dz;
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
                    CI(1:length(indices_m)) .* ( (Dxm+Dym+Dzm)/3 + Fractures(f1).Thickness ) ./...
                    ( ( (Dxm+Dym+Dzm)/3 ./ Reservoir.K(indices_m,1) ) + ( Fractures(f1).Thickness ./ Fractures(f1).K(g1,1) ) );
   
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
        function AddHarmonicConductivities(obj, Reservoir, Fractures)
            Dxm = obj.ReservoirGrid.dx;
            Dym = obj.ReservoirGrid.dy;
            Dzm = obj.ReservoirGrid.dz;
            Nm = obj.ReservoirGrid.N;
            for If1Local = 1 : length(obj.CrossConnections)
                Cells = obj.CrossConnections(If1Local).Cells;
                CI = obj.CrossConnections(If1Local).CI;
                If1Global = Nm+If1Local; % Global index of this fracture cell;
                Ind_frac1_Local = obj.Index_Global_to_Local(If1Global);
                f1 = Ind_frac1_Local.f;
                g1 = Ind_frac1_Local.g;
                
                indices_m = Cells( Cells <= Nm );
                obj.CrossConnections(If1Local).T_Geo_Cond(1:length(indices_m)) = ...
                    CI(1:length(indices_m)) .* ( (Dxm+Dym+Dzm)/3 + Fractures(f1).Thickness ) ./...
                    ( ( (Dxm+Dym+Dzm)/3 ./ Reservoir.K_Cond_eff(indices_m,1) ) + ( Fractures(f1).Thickness ./ Fractures(f1).K_Cond_eff(g1,1) ) );
   
                indices_f = Cells( Cells > Nm );
                if ~isempty(indices_f)
                    for n = 1:length(indices_f)
                        If2Global = indices_f(n);
                        If2Local = If2Global - Nm;
                        Ind_frac2_Local = obj.Index_Global_to_Local(indices_f(n));
                        f2 = Ind_frac2_Local.f;
                        g2 = Ind_frac2_Local.g;
                        obj.CrossConnections(If1Local).T_Geo_Cond(length(indices_m)+n) = ...
                            CI(length(indices_m)+n) * ( (obj.FracturesGrid.Grids(f1).dx + obj.FracturesGrid.Grids(f1).dy)/2 + Fractures(f2).Thickness ) ./...
                              ( ( (obj.FracturesGrid.Grids(f1).dx + obj.FracturesGrid.Grids(f1).dy)/2 ./ Fractures(f1).K_Cond_eff(g1,1) ) + ( Fractures(f2).Thickness ./ Fractures(f2).K_Cond_eff(g2,1) ) );
                    end
                end
            end
        end
        function AverageMassOnCoarseBlocks(obj, Status, FluidModel, Formulation)
            % virtual call
        end
    end
end