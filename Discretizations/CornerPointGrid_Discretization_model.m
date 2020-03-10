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
                PerfList = [];
                for p = 1:size(Wells.Inj(w).Coord)-1
                    PointA = Wells.Inj(w).Coord(p  ,:);
                    PointB = Wells.Inj(w).Coord(p+1,:);
                    PointM = (PointA+PointB)/2;
                    [ ~ , indList ] = min( vecnorm(PointM - obj.CornerPointGridData.Cell.Centroid,2,2 ) );
                    Count = 1;
                    while Count <= length(indList)
                        I = indList(Count);
                        Hexahedron.NW_Top_Corner = obj.CornerPointGridData.Cell(I).NW_Top_Corner;
                        Hexahedron.NE_Top_Corner = obj.CornerPointGridData.Cell(I).NE_Top_Corner;
                        Hexahedron.SW_Top_Corner = obj.CornerPointGridData.Cell(I).SW_Top_Corner;
                        Hexahedron.SE_Top_Corner = obj.CornerPointGridData.Cell(I).SE_Top_Corner;
                        Hexahedron.NW_Bot_Corner = obj.CornerPointGridData.Cell(I).NW_Bot_Corner;
                        Hexahedron.NE_Bot_Corner = obj.CornerPointGridData.Cell(I).NE_Bot_Corner;
                        Hexahedron.SW_Bot_Corner = obj.CornerPointGridData.Cell(I).SW_Bot_Corner;
                        Hexahedron.SE_Bot_Corner = obj.CornerPointGridData.Cell(I).SE_Bot_Corner;
                        Hexahedron.N_Faces = 6;
                        Hexahedron.Face(1).Pos = 'North';
                        Hexahedron.Face(1).PointA = Hexahedron.NW_Top_Corner;
                        Hexahedron.Face(1).PointB = Hexahedron.NW_Bot_Corner;
                        Hexahedron.Face(1).PointC = Hexahedron.NE_Bot_Corner;
                        Hexahedron.Face(1).PointD = Hexahedron.NE_Top_Corner;
                        AB = PointB - PointA;
                        PointC = obj.CornerPointGridData.Cell.Centroid;
                        AC = PointC - PointA;
                        Cos_Theta = dot( AB , L2B-L1A ) / ( norm(L1B-L1A) * norm(L2B-L1A) );
                    end
                end
            end
            
            % Producers
            for w = 1:Wells.NofProd
                I = Wells.Prod(w).Coord(1,1):1:Wells.Prod(w).Coord(2,1);
                J = Wells.Prod(w).Coord(1,2):1:Wells.Prod(w).Coord(2,2);
                K = Wells.Prod(w).Coord(1,3):1:Wells.Prod(w).Coord(2,3);
                if (sum(I > obj.ReservoirGrid.Nx) == 1 || sum(J > obj.ReservoirGrid.Ny) == 1 || sum(K > obj.ReservoirGrid.Nz) == 1)
                    error(['ERROR: coordinates of producer num ', num2str(w),' fall outside of the domain']);
                else
                     Wells.Prod(w).Cells = I + (J-1)*obj.ReservoirGrid.Nx + (K-1)*obj.ReservoirGrid.Nx*obj.ReservoirGrid.Ny;
                     Wells.Prod(w).ResizeObjects(length(Wells.Prod(w).Cells));
                end
            end
        end
    end
end