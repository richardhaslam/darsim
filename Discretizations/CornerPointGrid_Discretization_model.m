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
                if ~isempty(indices_m)
                    obj.CrossConnections(If1Local).T_Geo(1:length(indices_m)) = ...
                        CI(1:length(indices_m)) .* ( obj.ReservoirGrid.Volume(indices_m).^(1/3) + Fractures(f1).Thickness ) ./...
                        ( ( obj.ReservoirGrid.Volume(indices_m).^(1/3) ./ Reservoir.K(indices_m,1) ) + ( Fractures(f1).Thickness ./ Fractures(f1).K(g1,1) ) );
                end
                
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
            Nm = obj.ReservoirGrid.N;
            for If1Local = 1 : length(obj.CrossConnections)
                Cells = obj.CrossConnections(If1Local).Cells;
                CI = obj.CrossConnections(If1Local).CI;
                If1Global = Nm+If1Local; % Global index of this fracture cell;
                Ind_frac1_Local = obj.Index_Global_to_Local(If1Global);
                f1 = Ind_frac1_Local.f;
                g1 = Ind_frac1_Local.g;
                
                indices_m = Cells( Cells <= Nm );
                if ~isempty(indices_m)
                    obj.CrossConnections(If1Local).T_Geo_Cond(1:length(indices_m)) = ...
                        CI(1:length(indices_m)) .* ( obj.ReservoirGrid.Volume(indices_m).^(1/3) + Fractures(f1).Thickness ) ./...
                        ( ( obj.ReservoirGrid.Volume(indices_m).^(1/3) ./ Reservoir.K_Cond_eff(indices_m,1) ) + ( Fractures(f1).Thickness ./ Fractures(f1).K_Cond_eff(g1,1) ) );
                end
                
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
                PointA = Well.Coordinate.Value(p  ,:);
                PointB = Well.Coordinate.Value(p+1,:);
                LineSegment = lineSegment_DARSim(PointA,PointB);
                
                [ ~ , indList ] = min( vecnorm(LineSegment.PointM - obj.CornerPointGridData.Cells.Centroid, 2,2) );
                Count = 1;
                while Count <= length(indList)
                    I = indList(Count);
                    Count = Count+1;
                    NW_Top = obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(I,1) , : );
                    SW_Top = obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(I,2) , : );
                    SE_Top = obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(I,3) , : );
                    NE_Top = obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(I,4) , : );
                    NW_Bot = obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(I,5) , : );
                    SW_Bot = obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(I,6) , : );
                    SE_Bot = obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(I,7) , : );
                    NE_Bot = obj.CornerPointGridData.Nodes( obj.CornerPointGridData.Cells.Vertices(I,8) , : );
                    
                    if isempty(obj.CornerPointGridData.Cells.Faces{I})
                        ReservoirCell = hexahedron_DARSim(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot);
                    else
                        for n = 1 : length(obj.CornerPointGridData.Cells.Faces{I})
                            InterfaceFullIndex = obj.CornerPointGridData.Cells.Faces{I}(n);
                            InterfaceIndex = find( obj.CornerPointGridData.Internal_Faces.FullIndex == InterfaceFullIndex );
                            if ~isempty(InterfaceIndex) % This is an internal face
                                NodesIndex = unique( obj.CornerPointGridData.Internal_Faces.Vertices{InterfaceIndex} , 'stable' );
                            else                        % This is an external face
                                InterfaceIndex = find( obj.CornerPointGridData.External_Faces.FullIndex == InterfaceFullIndex );
                                NodesIndex = unique( obj.CornerPointGridData.External_Faces.Vertices{InterfaceIndex} , 'stable' );
                            end
                            Vertices = obj.CornerPointGridData.Nodes( NodesIndex , : );
                            n_vec = cross( Vertices(2,:)-Vertices(1,:) , Vertices(end,:)-Vertices(1,:) );
                            Face(n,1) = polygon_DARSim( length(NodesIndex) , n_vec , mean(Vertices) );
                            Face(n,1).AddVertices(Vertices);
                            Face(n,1).Centroid = mean(Vertices);
                        end
                        ReservoirCell = hexahedron_MRST(NW_Top,SW_Top,SE_Top,NE_Top,NW_Bot,SW_Bot,SE_Bot,NE_Bot);
                        ReservoirCell.AddFaces(Face);
                    end
                    ReservoirCell.Centroid = obj.CornerPointGridData.Cells.Centroid(I,:);
                    Epsilon = 1e-10 * ( min(obj.ReservoirGrid.Volume) )^(1/3);
                    [Geostatus, IntersectPoints] = ReservoirCell.Obtain_Polyhedron_LineSegment_Intersection(LineSegment, Epsilon);
                    
                    % If it is the first try and no intersection occurs between the well and
                    % this reservoir cell, add the neighboring cells to the list for intersection check.
                    if isempty(IntersectPoints)
                        if Count == 2
                            indNeighbors = obj.CornerPointGridData.Cells.Neighbors{I};
                            indList = union(indList, indNeighbors, 'stable');
                        end
                        continue;
                    else
                        % There was already an intersection between the well and this reservoir cell.
                        % Therefore, add the neighboring cells to the list for intersection check
                        % as the neighbors may intersect as well.
                        indNeighbors = obj.CornerPointGridData.Cells.Neighbors{I};
                        indList = union(indList, indNeighbors, 'stable');
                    end
                    
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