% class plate for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2018-04-20
%Last modified: 2018-04-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef plate_FracGen < handle
    properties
        PointA
        PointB
        PointC
        PointD
        PointM
        Points
        Area
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
    end
end