% Thermal compressible phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo ...
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Geothermal_SinglePhase_phase < phase
    properties
        rho0 % Reference density
        mu0 % Reference viscosity
        cf0 % Reference fluid compressibility
        Kf % Conductivity
        Cp_std % Specific Heat
        uws = 420000;% internal energy at saturation J/kg
        Tsat = 373; % T at saturation condition (assumed constant 100 C)
        Psat = 1e5; % P at saturation condition (assumed 1e5 Pa)
    end
    methods
        function cond = AddConductivity(obj)
            cond = obj.Kf;
        end
        function rho = ComputeDensityBasedOnTemperature(obj, p, T)
            cw =  (0.0839*T.^2 - 64.593*T + 12437)*1e-12;
            rhofs = -0.0032*T.^2 + 1.7508*T + 757.5; 
            rho = rhofs.*(1+cw.*(p-obj.Psat));
        end
        function rho = ComputeDensityBasedOnEnthalpy(obj, i, PhaseIndex, p, h)
            % The index "i" is "1" for water and "2" for steam.
            for k = 1:3
                if i == 1                    
                    rho(PhaseIndex == k) = ( ...
                        1.00207 + ...
                        4.42607e-11.*(p(PhaseIndex == k).*1e1) + ...
                        -5.47456e-12.*(h(PhaseIndex == k).*1e4) + ...
                        5.02875e-21.*(h(PhaseIndex == k).*1e4).*(p(PhaseIndex == k).*1e1) + ...
                        -1.24791e-21.*(h(PhaseIndex == k).*1e4).^2 ...
                        ).*1e3;
                elseif i == 2                    
                    rho(PhaseIndex == k) = ( ...
                        -2.26162e-5 + ...
                        4.38441e-9.*(p(PhaseIndex == k).*1e1) + ...
                        -1.79088e-19.*(p(PhaseIndex == k).*1e1).*(h(PhaseIndex == k).*1e4) + ...
                        3.69276e-36.*(p(PhaseIndex == k).*1e1).^4 + ...
                        5.17644e-41.*(p(PhaseIndex == k).*1e1).*(h(PhaseIndex == k).*1e4).^3 ...
                        ).*1e3;
                end
                rho = rho';
            end
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
            %
            B = 247.8;
            %
            C = T-140;
            dCdT = 1;
            d2CdT2 = 0;
            %
            D = B./C;
            dDdC = -B./C.^2;
            d2DdC2 = (2.*B./C.^3);
            %
            E = 10.^D;
            dEdD = E .* log(D);
            d2EdD2 = dEdD .* log(D) + E ./ D;
            %
            mu = A.*E;
            dmudE = A;
            d2mudE2 = 0;
            %
            dmudT = dmudE.*dEdD.*dDdC.*dCdT;
            d2mudT2 = d2mudE2.*dEdD.*dDdC.*dCdT + dmudE.*d2EdD2.*dDdC.*dCdT + dmudE.*dEdD.*d2DdC2.*dCdT + dmudE.*dEdD.*dDdC.*d2CdT2;
        end
        function h = ComputeEnthalpy(obj, p, T)
            rho = obj.ComputeDensityBasedOnTemperature(p, T);
            h = obj.uws + obj.Cp_std*(T-obj.Tsat)+p./rho;
        end
        function dhdp = ComputeDhDp(obj, p, T)
            rho = obj.ComputeDensityBasedOnTemperature(p, T);
            drhodp = obj.ComputeDrhoDp(p, T);
            dhdp = ((rho - p.*drhodp)./rho.^2);
        end
        function [dhdT,d2hdT2] = ComputeDhDT(obj, p, T)
            rho = obj.ComputeDensityBasedOnTemperature(p, T);
            [drhodT,d2rhodT2] = obj.ComputeDrhoDT(p, T);
            dhdT = obj.Cp_std + p.*(-drhodT./rho.^2);
            d2hdT2 = p.* ( -d2rhodT2./rho.^2 + 2.*drhodT./rho.^3 );
        end
        % Injection properties; we are injecting only water, so it is
        % easier to use existing functions from Geothermal singlephase
        function rho = ComputeWaterDensity(obj, p, T)
            cw =  (0.0839*T.^2 - 64.593*T + 12437)*1e-12;
            rhofs = -0.0032*T.^2 + 1.7508*T + 757.5; 
            rho = rhofs.*(1+cw.*(p-obj.Psat));
        end
        function h = ComputeWaterEnthalpy(obj, p, T)
            A = -2.41231;
            B = 2.5622e-8;
            C = -9.31415e-17;
            D = -2.2568e-19;
            % This is the re-ordered formula from the "ComputeTemperature" function above in the class "Geothermal_Multiphase_phase".
            h = ( - B + sqrt( B^2 - 4*D*(A+C*(p.*1e1).^2-(T-273.15)) ) ) / (2*D*1e4);
        end
        function T = ComputeWaterTemperature(obj, p, h)
            i = 1; % the "i=1" is for water phase
            PhaseIndex = 1; % the phase index "1" refers to water phase
            rho = obj.ComputeDensityBasedOnEnthalpy(i, PhaseIndex, p, h);
            T = ( (h - obj.uws - p./rho) / obj.Cp_std ) + obj.Tsat;
        end
        function mu = ComputeWaterViscosity(obj, T)
            A = 2.414e-5;   B = 247.8;  C = T-140;   D = B./C;   E = 10.^D;            
            mu = A.*E;
        end
        function v = ComputeVelocity(obj, p, mu)
            % virtual call
        end
    end
end