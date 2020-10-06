% Thermal compressible phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Arjan Marelis
%TU Delft
%Created: 21 January 2020
%Last modified: 23 January 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Geothermal_Multiphase_phase < phase
    properties
        Kf
        Cp_std              % Specific Heat of Phase in standard condition
        uws = 420000;       % internal energy at saturation J/kg
        Tsat = 373;         % T at saturation condition (assumed constant 100 C)
        Psat = 1e5;         % P at saturation condition (assumed 1e5 Pa)
    end
    methods
        function cond = AddConductivity(obj)
            cond = obj.Kf;
        end
        function PhaseEnthalpy = ComputePhaseEnthalpy(obj, i, p)
            if i == 1
                PhaseEnthalpy = (7.30984e9 + 1.29239e2.*(p.*1e1) - 1.00333e-6.*(p.*1e1).^2 + 3.9881e-15.*(p.*1e1).^3 + ...
                                 - 9.90697e15.*(p.*1e1).^-1 + 1.29267e22.*(p.*1e1).^-2 - 6.28359e27.*(p.*1e1).^-3).*1e-7 .*1e3;
            elseif i == 2
                PhaseEnthalpy = (2.82282e10 - 3.91952e5.*(p.*1e1).^-1 + 2.54342e21.*(p.*1e1).^-2 - 9.38879e-8.*(p.*1e1).^2).*1e-7 .*1e3;
            end
        end
        function T = ComputeTemperature(obj, PhaseIndex, p, h)
            % Compressed water region [K]
            T(PhaseIndex == 1) = 273.15 - 2.41231 + 2.56222e-8.*(h(PhaseIndex == 1).*1e4) - 9.31415e-17.*(p(PhaseIndex == 1).*1e1).^2 - 2.2568e-19.*(h(PhaseIndex == 1).*1e4).^2;
            % Twophase region
            T(PhaseIndex == 2) = 273.15 - 2.41231 + 2.56222e-8.*(h(PhaseIndex == 2).*1e4) - 9.31415e-17.*(p(PhaseIndex == 2).*1e1).^2 - 2.2568e-19.*(h(PhaseIndex == 2).*1e4).^2;
            % Superheated steam region
            T(PhaseIndex == 3) = 273.15 - 374.669 + 4.79921e-6.*(p(PhaseIndex == 3).*1e1) - 6.33606e-15.*(p(PhaseIndex == 3).*1e1).^2 + ...
                7.39386e-19.*(h(PhaseIndex == 3).*1e4).^2 - 3.3372e34.*(p(PhaseIndex == 3).*1e1).^-2.*(h(PhaseIndex == 3).*1e4).^-2 + ...
                3.57154e19.*(p(PhaseIndex == 3).*1e1).^-3 - 1.1725e-37.*(p(PhaseIndex == 3).*1e1).*(h(PhaseIndex == 3).*1e4).^3 + ...
                -2.26861e43.*(h(PhaseIndex == 3).*1e4).^-4;
            T = T';
        end
        function mu = ComputeViscosity(obj, i, PhaseIndex, T)
            for k = 1:3
                if i == 1
                    mu(PhaseIndex == k) = (241.4 .* 10.^(247.8./((T(PhaseIndex == k)-273.15) + 133.15)) ) .* 1e-4 .* 1e-3;
                elseif i == 2
                    mu(PhaseIndex == k) = (0.407.*(T(PhaseIndex == k)-273.15) + 80.4) .* 1e-4 .* 1e-3;
                end
            end
            mu = mu'; 
        end
        function psat = ComputeSaturationPressure(obj, Status, PhaseIndex)
            T = Status.Properties('T').Value;
            % saturation pressure
            theta(PhaseIndex == 2) = T(PhaseIndex == 2) + (-0.23855557567849 ./ (T(PhaseIndex == 2) - 0.65017534844798e3));
            A(PhaseIndex == 2) = theta(PhaseIndex == 2).^2 + 0.11670521452767e4.*theta(PhaseIndex == 2) + -0.72421316703206e6;
            B(PhaseIndex == 2) = -0.17073846940092e2.*theta(PhaseIndex == 2).^2 + 0.12020824702470e5.*theta(PhaseIndex == 2) + -0.32325550322333e7;
            C(PhaseIndex == 2) = 0.14915108613530e2.*theta(PhaseIndex == 2).^2 + -0.48232657361591e4.*theta(PhaseIndex == 2) + 0.40511340542057e6;
            % P = Psat [Pa]
            psat(PhaseIndex == 2) = ( (2.*C(PhaseIndex == 2)) ./ (-1.*B(PhaseIndex == 2)+sqrt(B(PhaseIndex == 2).^2-4.*A(PhaseIndex == 2).*C(PhaseIndex == 2))) ).^4 .* 1e6;
            % apply correction to size of psat for logical indexing of pressure values
            psat(psat == 0) = [];
        end
        function rho = ComputeDensity(obj, i, PhaseIndex, p, h)
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
        %% Derivatives (directly from tables)
        function drhodp = ComputeDrhoDp(obj, i, PhaseIndex, p, h)
            % Note that the derivative for two-phase region is implemented identical to single -phase regions; 
            % when Psat has no derivative, the equation changes and this function should change as well !
            for k = 1:3
                if i == 1                    
                    drhodp(PhaseIndex == k) = ( ...
                        4.42607e-11 + ...
                        5.02875e-21.*(h(PhaseIndex == k).*1e4) ...
                        ).*1e3;
                elseif i == 2                    
                    drhodp(PhaseIndex == k) = ( ...
                        4.38441e-9 + ...
                        - 1.79088e-19.*(h(PhaseIndex == k).*1e4) + ...
                        3.69276e-36.*4.*(p(PhaseIndex == k).*1e1).^3 + ...
                        5.17644e-41.*(h(PhaseIndex == k).*1e4).^3 ...
                        ).*1e3;
                end
                %drhodp = drhodp';
            end
        end
        function drhodh = ComputeDrhoDh(obj, i, PhaseIndex, p, h)
            % Note that the derivative for two-phase region is implemented identical to single -phase regions; 
            % when Psat has no derivative, the equation changes and this function should change as well !
            for k = 1:3
                if i == 1                    
                    drhodh(PhaseIndex == k) = ( ...
                        -5.47456e-12 + ...
                        5.02875e-21.*(p(PhaseIndex == k).*1e1) + ...
                        - 1.24791e-21.*2.*(h(PhaseIndex == k).*1e4) ...
                        ).*1e3;
                elseif i == 2                    
                    drhodh(PhaseIndex == k) = ( ...
                        -1.79088e-19.*(p(PhaseIndex == k).*1e1) + ...
                        5.17644e-41.*(p(PhaseIndex == k).*1e1).*3.*(h(PhaseIndex == k).*1e4).^2 ...
                        ).*1e3;
                end
                %drhodh = drhodh';
            end
        end
        function dTdp = ComputeDTDp(obj, PhaseIndex, p, h) 
            % Compressed water region [K]
            dTdp(PhaseIndex == 1) = -9.31415e-17.*2.*(p(PhaseIndex == 1).*1e1);
            % Twophase region
            dTdp(PhaseIndex == 2) = -9.31415e-17.*2.*(p(PhaseIndex == 2).*1e1);
            % Superheated steam region
            dTdp(PhaseIndex == 3) = 4.79921e-6 - 6.33606e-15.*2.*(p(PhaseIndex == 3).*1e1) + ...
                                     3.3372e34.*(h(PhaseIndex == 3).*1e4).^-2.*(p(PhaseIndex == 3).*1e1).^-3 + ...
                                     -3.57154e19.*3.*(p(PhaseIndex == 3).*1e1).^-4 - 1.1725e-37.*(h(PhaseIndex == 3).*1e4).^3;
            dTdp = dTdp';
        end
        function dTdh = ComputeDTDh(obj, PhaseIndex, p, h)
            % Compressed water region [K]
            dTdh(PhaseIndex == 1) = 2.56222e-8 - 2.2568e-19.*2.*(h(PhaseIndex == 1).*1e4);
            % Twophase region
            dTdh(PhaseIndex == 2) = 2.56222e-8 - 2.2568e-19.*2.*(h(PhaseIndex == 2).*1e4);
            % Superheated steam region
            dTdh(PhaseIndex == 3) = 7.39386e-19.*2.*(h(PhaseIndex == 3).*1e4) + ...
                                     3.3372e34.*2.*(h(PhaseIndex == 3).*1e4).^-3.*(p(PhaseIndex == 3).*1e1).^-2 + ...
                                     -1.1725e-37.*3.*(h(PhaseIndex == 3).*1e4).^2.*(p(PhaseIndex == 3).*1e1) + ...
                                     2.26861e43.*4.*(h(PhaseIndex == 3).*1e4).^-5;
            dTdh = dTdh';
        end
        % These depend on how we treat the conductive flux term
        function d2Td2p = ComputeD2TD2p(obj, Pgrid, Hgrid, TTable, h, p) 
            %
        end
        function d2Td2h = ComputeD2TD2h(obj, Pgrid, Hgrid, TTable, h, p)
            %
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
            % This is the re-ordered formula from the "ComputeTemperature" function above in this class.
            h = ( - B + sqrt( B^2 - 4*D*(A+C*(p.*1e1).^2-(T-273.15)) ) ) / (2*D*1e4);
        end
        function T = ComputeWaterTemperature(obj, p, h)
            PhaseIndex = 1; % the phase index "1" refers to water phase
            T = obj.ComputeTemperature(PhaseIndex, p, h);
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
