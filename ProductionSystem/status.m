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
    p
    pc
    S
    x
    z
    rho
    rhoT
    ni
end
methods
    function Initialize(obj, N, nf, ncomp)
        % Create objects
        obj.p = ones(N, 1);
        obj.pc = zeros(N, 1);
        obj.z = ones(N, ncomp);
        obj.S = zeros(N, 1);
        obj.ni = 0.5*ones(N, 1);
        obj.x = zeros(N, nf*ncomp);
        obj.rho = zeros(N, nf);
        obj.rhoT = zeros(N, 1); 
    end
end
end