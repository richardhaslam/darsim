% Fracture system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fractures_mousa < handle
    properties
        Length_AB
        Width_AD
        Aperture
        Porosity
        Permeability
        PointA
        PointB
        PointC
        PointD
        PointM
        Points
        Area_Total
        AB_vec
        AD_vec
        AC_vec
        BC_vec
        CD_vec
        DA_vec
        n_vec
        n_vec_plot
        Equation
        isParallel
        N_Length_AB
        N_Width_AD
        ADM
        D_Length_AB
        D_Width_AD
        D_Length_AB_vec
        D_Width_AD_vec
        CellCenterCoords
        CellCenterCoordsV2
        GridCoords
        intersectCoord_matCell
        areaFrac_matCell
        areaFrac_matCell_sum
        aveDist_matCell
        T_Geo_matCell
        intersectCoord_fracObj
        intersectCoord_fracCell
        Overlap_frac
        areaFrac_fracCell
        areaFrac_fracCell_sum
        aveDist_fracCell
        T_Geo_fracCell
        NumOf_fracCellConn
        StarDelta2D_Neighbor
    end
    
    
end