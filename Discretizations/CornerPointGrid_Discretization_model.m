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
        function DefinePerforatedCells(obj, Wells)
            % Injectors
            for w = 1:Wells.NofInj
                for p = 1:size(Wells.Inj(w).Coord,1) - 1
                    PointA = Wells.Inj(w).Coord(p  ,:)';
                    PointB = Wells.Inj(w).Coord(p+1,:)';
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
                            Wells.Inj(w).Cells = [ Wells.Inj(w).Cells ; I ];
                        end
>>>>>>> ee71a6b9ad3950abcf7ce630eebf8898126a1d01
                    end
                end
                Wells.Inj(w).Cells  = sort(Wells.Inj(w).Cells);
                Wells.Inj(w).ResizeObjects(length(Wells.Inj(w).Cells));
            end
            
            % Producers
            for w = 1:Wells.NofProd
                for p = 1:size(Wells.Prod(w).Coord,1) - 1
                    PointA = Wells.Prod(w).Coord(p  ,:)';
                    PointB = Wells.Prod(w).Coord(p+1,:)';
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
                            Wells.Prod(w).Cells = [ Wells.Prod(w).Cells ; I ];
                        end
                    end
                end
                Wells.Prod(w).Cells  = sort(Wells.Prod(w).Cells);
                Wells.Prod(w).ResizeObjects(length(Wells.Prod(w).Cells));
            end
        end
    end
end