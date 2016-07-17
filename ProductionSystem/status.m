% Status 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 17 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef status < handle
properties
    Pot
    p
    S
    x1
    z
    rho
    rhoT
end
methods
    function Initialize(obj, N)
        % Define initial values
        P_init = linspace(1e5, 10e4, N);
        z_init = ones(N, 1)*0.0;
        
        % Create objects
        obj.p = ones(N, 1);
        obj.z = ones(N, 2);
        obj.S = zeros(N, 1);
        obj.x1 = zeros(N, 2);
        obj.rho = zeros(N, 2);
        obj.rhoT = zeros(N, 1);
        
        % Assign initial valus
        obj.p = obj.p .* P_init';
        obj.z(:,1) = z_init;      
        obj.z(:,2) = 1 - z_init;   
    end
end
end