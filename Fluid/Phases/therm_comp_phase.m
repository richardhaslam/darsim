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
        function cond = AddConductivity(obj, p, T)
            cond = obj.Kf * ones(size(p));
        end
        function rho = ComputeDensity(obj, p, T)
            cw =  (0.0839*T.^2 - 64.593*T + 12437)*1e-12;
            rhofs = -0.0032*T.^2 + 1.7508*T + 757.5; 
            rho = rhofs.*(1+cw.*(p-obj.Psat));
        end
        function drhodp = ComputeDrhoDp(obj, p, T)
            cw =  (0.0839*T.^2 - 64.593*T + 12437)*1e-12;
            rhofs = -0.0032*T.^2 + 1.7508*T + 757.5;
            drhodp = rhofs.*cw;
        end
        function [drhodT,d2rhodT2] = ComputeDrhoDT(obj, p, T)
            cw =  (0.0839*T.^2 - 64.593*T + 12437)*1e-12;
            dcwdT = (2*0.0839*T - 64.593)*1e-12;
            d2cwdT2 = 2*0.0839*1e-12;
            rhofs = -0.0032*T.^2 + 1.7508*T + 757.5;
            drhofsdT = 2*(-0.0032).*T + 1.7508;
            d2rhofsdT2 = 2*(-0.0032);
            drhodT = drhofsdT + (p-obj.Psat).*(rhofs.*dcwdT + cw.*drhofsdT);
            d2rhodT2 = d2rhofsdT2 + (p-obj.Psat).* ( rhofs.*d2cwdT2 + drhofsdT.*dcwdT + cw.*d2rhofsdT2 + dcwdT.*drhofsdT);
        end
        function [mu,dmudT,d2mudT2] = ComputeViscosity(obj, T)
            A = 2.414e-5;
            
            B = 247.8;
            
            C = T-140;
            dCdT = 1;
            d2CdT2 = 0;
            
            D = B./C;
            dDdC = -B./C.^2;
            d2DdC2 = (2.*B./C.^3);
            
            E = 10.^D;
            dEdD = E .* log(D);
            d2EdD2 = dEdD .* log(D) + E ./ D;
            
            mu = A.*E;
            dmudE = A;
            d2mudE2 = 0;
            
            dmudT = dmudE.*dEdD.*dDdC.*dCdT;
            d2mudT2 = d2mudE2.*dEdD.*dDdC.*dCdT + dmudE.*d2EdD2.*dDdC.*dCdT + dmudE.*dEdD.*d2DdC2.*dCdT + dmudE.*dEdD.*dDdC.*d2CdT2;
        end
        function h = ComputeEnthalpy(obj, p, T)
            rho = obj.ComputeDensity(p, T);
            h = obj.uws + obj.Cp*(T-obj.Tsat)+p./rho;
        end
        function dhdp = ComputeDhDp(obj, p, T)
            rho = obj.ComputeDensity(p, T);
            drhodp = obj.ComputeDrhoDp(p, T);
            dhdp = ((rho - p.*drhodp)./rho.^2);
        end
        function [dhdT,d2hdT2] = ComputeDhDT(obj, p, T)
            rho = obj.ComputeDensity(p, T);
            [drhodT,d2rhodT2] = obj.ComputeDrhoDT(p, T);
            dhdT = obj.Cp + p.*(-drhodT./rho.^2);
            d2hdT2 = p.* ( -d2rhodT2./rho.^2 + 2.*drhodT./rho.^3 );
        end
        function v = ComputeVelocity(obj, p, mu)
            % virtual call
        end
    end
end