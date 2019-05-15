% Thermal compressible phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo ...
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef therm_comp_phase < phase
    properties
        rho0 % Reference density
        mu0 % Reference viscosity
        cf0 % Reference fluid compressibility
        Kf % Conductivity
        Cp % Specific Heat
        uws = 420000;% internal energy at saturation J/kg
        Tsat = 373; % T at saturation condition (assumed constant 100 C)
        Psat = 1e5; % P at saturation condition (assumed 1e5 Pa)
    end
    methods
        function rho = ComputeDensity(obj, p, T)
            cf =  (0.0839*T.^2 - 64.593*T + 12437)*1e-12;
            rho_fs = -0.0032*T.^2 + 1.7508*T + 757.5; 
            rho = rho_fs.*(1+cf.*(p-obj.Psat));
        end
        function drhodp = ComputeDrhoDp(obj, p, T)
            cf =  (0.0839*T.^2 - 64.593*T + 12437)*1e-12;
            rho_fs = -0.0032*T.^2 + 1.7508*T + 757.5;
            drhodp = rho_fs.*cf;
        end
        function drhodT = ComputeDrhoDT(obj, p, T)
            cf =  (0.0839*T.^2 - 64.593*T + 12437)*1e-12;
            dcwdT = (2*0.0839*T - 64.593)*1e-12;
            rho_fs = -0.0032*T.^2 + 1.7508*T + 757.5;
            drho_fs_dT = 2*-0.0032.*T + 1.7508;
            drhodT = drho_fs_dT + (p-obj.Psat).*(rho_fs.*dcwdT + cf.*drho_fs_dT);
        end
        function [mu,dmudT] = ComputeViscosity(obj, T)
%             mu = 1.5396e-3*exp(-0.018.*(T-273));
%             dmudT = -0.018*1.5396e-3*exp(-0.018.*(T-273));

            A = 2.414e-5;
            B = 247.8;
            C = T-140; dCdT = 1;
            D = B./C; dDdC = -B./C.^2;
            E = 10.^D; dEdD = 10.^D.*log(D);
            mu = A.*E; dmudE = A;           
            dmudT = 0*dmudE.*dEdD.*dDdC.*dCdT;
              
%               A = 0.01779;
%               B = (1 + 0.03368.*T + 0.000220099.*T.^2); 
%               dBdT = 0.03368 + 0.000220099.*T;
%               mu = A./B;
%               dmudT = -A./B.^2 .* dBdT;
        end
%         function dmudT = ComputeDmu(obj, T)
%             [~,dmudT] = ComputeViscosity(obj, T);
%         end
        function h = ComputeEnthalpy(obj, p, T)
            rho = obj.ComputeDensity(p, T);
            h = obj.uws + obj.Cp*(T-obj.Tsat)+p./rho;
        end
        function v = ComputeVelocity(obj, p, mu)
%             virtual call
        end
        function dhdp = ComputeDhDp(obj, p, T)
            rho = obj.ComputeDensity(p, T);
            drhodp = obj.ComputeDrhoDp(p, T);
            dhdp = ((rho - p.*drhodp)./rho.^2);
        end
        function dhdT = ComputeDhDT(obj, p, T)
            rho = obj.ComputeDensity(p, T);
            drhodT = obj.ComputeDrhoDT(p, T);
            dhdT = obj.Cp + p.*(-drhodT./rho.^2);
        end
    end
end