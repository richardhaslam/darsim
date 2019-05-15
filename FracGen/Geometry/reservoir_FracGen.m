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
        alpha_Tx
        alpha_Ty
        alpha_Tz
        % alpha is a the fraction of blockage (from 0 to 1) between
        % neighboring matrix cells due to existence of fracture planes when
        % considering pEDFM method. The alpha=0 means no blockage and
        % alpha=1 means total blockage.
        Epsilon
        OvelapList = [];
    end
    methods
        function Discretize(obj, CompareAccuracy, Type)
            obj.DX  = obj.LX/obj.NX;
            obj.DY  = obj.LY/obj.NY;
            obj.DZ  = obj.LZ/obj.NZ;
            obj.Xcm = linspace(obj.DX/2, obj.LX - obj.DX/2, obj.NX)';              % Grid centers locations in x-direction
            obj.Ycm = linspace(obj.DY/2, obj.LY - obj.DY/2, obj.NY)';              % Grid centers locations in y-direction
            obj.Zcm = linspace(obj.DZ/2, obj.LZ - obj.DZ/2, obj.NZ)';              % Grid centers locations in z-direction
            obj.Xim = linspace(0, obj.LX, obj.NX+1)';                               % Interface locations in x-direction
            obj.Yim = linspace(0, obj.LY, obj.NY+1)';                               % Interface locations in y-direction
            obj.Zim = linspace(0, obj.LZ, obj.NZ+1)';                               % Interface locations in z-direction
            obj.Epsilon = CompareAccuracy * min([obj.DX , obj.DY , obj.DZ]);
            if strcmp(Type,'pEDFM')
                obj.alpha_Tx = zeros(obj.NX+1, obj.NY  , obj.NZ  );
                obj.alpha_Ty = zeros(obj.NX  , obj.NY+1, obj.NZ  );
                obj.alpha_Tz = zeros(obj.NX  , obj.NY  , obj.NZ+1);
            end
        end
        function PrintInfo(obj, Domain)
            disp( char(5));
            if strcmp(Domain,'2D')
                disp( '2D Reservoir geometry:');
                disp( ['Length: ', num2str(obj.LX), ' [m]'] );
                disp( ['Width : ', num2str(obj.LY), ' [m]'] );
                disp( ['Grid  : ', num2str(obj.NX), ' x ',  num2str(obj.NY),...
                      ' = ', num2str(obj.NX * obj.NY)] );
            else
                disp( '3D Reservoir geometry:');
                disp( ['Length: ', num2str(obj.LX), ' [m]'] );
                disp( ['Width : ', num2str(obj.LY), ' [m]'] );
                disp( ['Depth : ', num2str(obj.LZ), ' [m]'] );
                disp( ['Grid  : ', num2str(obj.NX), ' x ',  num2str(obj.NY), ' x ', num2str(obj.NZ),...
                      ' = ', num2str(obj.NX * obj.NY * obj.NZ)] );
            end
            disp( '---------------------------------------------------------' );
        end
    end
end