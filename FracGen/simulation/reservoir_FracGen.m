% class reservoir for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini & Mousa HosseiniMehr
%TU Delft
%Created: 2016-07-04
%Last modified: 2017-03-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reservoir_FracGen < handle
    properties
        LX
        LY
        LZ
        NX
        NY
        NZ
        DX
        DY
        DZ
        Xcm
        Ycm
        Zcm
        Xim
        Yim
        Zim
        almostZero
    end
    methods
        function Discretize(obj, CompareAccuracy, Dimension_Type)
            if Dimension_Type == '3D'
                obj.DX  = obj.LX/obj.NX;
                obj.DY  = obj.LY/obj.NY;
                obj.DZ  = obj.LZ/obj.NZ;
                obj.Xcm = linspace(obj.DX/2, obj.LX - obj.DX/2, obj.NX);               % Grid centers locations in x-direction
                obj.Ycm = linspace(obj.DY/2, obj.LY - obj.DY/2, obj.NY);               % Grid centers locations in y-direction
                obj.Zcm = linspace(obj.DZ/2, obj.LZ - obj.DZ/2, obj.NZ);               % Grid centers locations in z-direction
                obj.Xim = linspace(0, obj.LX, obj.NX+1);                               % Interface locations in x-direction
                obj.Yim = linspace(0, obj.LY, obj.NY+1);                               % Interface locations in y-direction
                obj.Zim = linspace(0, obj.LZ, obj.NZ+1);                               % Interface locations in z-direction
                obj.almostZero = CompareAccuracy * min([obj.DX , obj.DY , obj.DZ]);
            elseif Dimension_Type == '2D'
                obj.DX  = obj.LX/obj.NX;
                obj.DY  = obj.LY/obj.NY;
                obj.DZ  = NaN;
                obj.Xcm = linspace(obj.DX/2, obj.LX - obj.DX/2, obj.NX);               % Grid centers locations in x-direction
                obj.Ycm = linspace(obj.DY/2, obj.LY - obj.DY/2, obj.NY);               % Grid centers locations in y-direction
                obj.Zcm = NaN;                                                         % Grid centers locations in z-direction
                obj.Xim = linspace(0, obj.LX, obj.NX+1);                               % Interface locations in x-direction
                obj.Yim = linspace(0, obj.LY, obj.NY+1);                               % Interface locations in y-direction
                obj.Zim = NaN;                                                         % Interface locations in z-direction
                obj.almostZero = CompareAccuracy * min([obj.DX , obj.DY]);
            end
        end
    end
end