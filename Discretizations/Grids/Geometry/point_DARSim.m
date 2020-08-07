% class point for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2019-02-07
%Last modified: 2019-02-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For now, this class will not be used.
classdef point_DARSim < handle
    properties
        X
        Y
        Z
    end
    methods
        function obj = point_DARSim(x,y,z)
            obj.X = x;
            obj.Y = y;
            obj.Z = z;
        end
        function newPoint = plus(obj,Point)
            x = obj.X + Point.X;
            y = obj.Y + Point.Y;
            z = obj.Z + Point.Z;
            newPoint = point_DARSim(x,y,z);
        end
        function newPoint = minus(obj,Point)
            x = obj.X - Point.X;
            y = obj.Y - Point.Y;
            z = obj.Z - Point.Z;
            newPoint = point_DARSim(x,y,z);
        end
        function PointAfterRotation = RotatePointAroundLine(obj, rotationAxis , pointInAxis , rotationAngle)
            rotationAxis = rotationAxis / norm( rotationAxis );
            Theta   = rotationAngle;
            
            a = pointInAxis(1); u = rotationAxis(1);
            b = pointInAxis(2); v = rotationAxis(2);
            c = pointInAxis(3); w = rotationAxis(3);
            
            rotationMatrix = [  u^2 + (v^2+w^2)*cos(Theta)         ,  u*v*(1-cos(Theta)) - w*sin(Theta)  ,  u*w*(1-cos(Theta)) + v*sin(Theta)  ,  ( a*(v^2+w^2) - u*(b*v+c*w) ) * (1-cos(Theta)) + (b*w-c*v)*sin(Theta)
                u*v*(1-cos(Theta)) + w*sin(Theta)  ,  v^2 + (u^2+w^2)*cos(Theta)         ,  v*w*(1-cos(Theta)) - u*sin(Theta)  ,  ( b*(u^2+w^2) - v*(a*u+c*w) ) * (1-cos(Theta)) + (c*u-a*w)*sin(Theta)
                u*w*(1-cos(Theta)) - v*sin(Theta)  ,  v*w*(1-cos(Theta)) + u*sin(Theta)  ,  w^2 + (u^2+v^2)*cos(Theta)         ,  ( c*(u^2+v^2) - w*(a*u+b*v) ) * (1-cos(Theta)) + (a*v-b*u)*sin(Theta)
                0                                  ,  0                                  ,   0                                 ,  1                                                                     ];
            
            ANS = rotationMatrix * [ obj.X ; obj.Y ; obj.Z ; 1 ];
            
            PointAfterRotation = point_FracGen( ANS(1) , ANS(2) , ANS(3) );
        end
    end
end