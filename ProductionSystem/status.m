% Status 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 17 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef status < matlab.mixin.Copyable
properties
    T
    Pot
    p
    pc
    S
    x1
    z
    rho
    rhoT
end
methods
    function Initialize(obj, N)
        % Create objects
        obj.p = zeros(N, 1);
        obj.pc = zeros(N, 1);
        obj.z = ones(N, 2);
        obj.S = zeros(N, 1);
        obj.x1 = zeros(N, 2);
        obj.rho = zeros(N, 2);
        obj.rhoT = zeros(N, 1); 
    end
end
end