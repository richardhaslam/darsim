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
                    
                    [Geostatus, IntersectPoints] = Hexahedron.Obtain_Hexahedron_LineSegment_Intersection(LineSegment);
                    
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