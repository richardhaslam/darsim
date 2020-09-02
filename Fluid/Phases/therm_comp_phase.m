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
        Pstepsize = 1e5;
        Tstepsize = 1; % implement this

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

        
%         function rho = ComputeDensity(obj, Pgrid, Tgrid, rhoTable, T, p)
%             rho = interp2(Tgrid, Pgrid, rhoTable, T, p, 'linear');
%         end
%         function mu = ComputeViscosity(obj, Pgrid, Tgrid, muTable, T, p)
%             mu = interp2(Tgrid, Pgrid, muTable, T, p, 'linear');
%         end
% %         function cond = AddConductivity(obj, Pgrid, Tgrid, ThermCondTable, T, p)
% %             cond = interp2(Tgrid, Pgrid, ThermCondTable, T, p, 'linear');
% %         end
%         function h = ComputeEnthalpy(obj, Pgrid, Tgrid, PhaseEnthalpyTable, T, p)
%             h = interp2(Tgrid, Pgrid, PhaseEnthalpyTable, T, p);
%         end
%         
%         function [dmudT,d2mudT2] = ComputeDmuDT(obj, Pgrid, Tgrid, muTable, T, p) 
%             [table_dmudT,~] = gradient(muTable,obj.Tstepsize,obj.Pstepsize);
%             dmudT = interp2(Tgrid, Pgrid, table_dmudT, T, p, 'linear');
%             [table_d2mudT2,~] = gradient(table_dmudT,obj.Tstepsize,obj.Pstepsize); 
%             d2mudT2 = interp2(Tgrid, Pgrid, table_d2mudT2, T, p, 'linear');   
%         end
%         function drhodp = ComputeDrhoDp(obj, Pgrid, Tgrid, rhoTable, T, p)
%             [~,table_drhodp] = gradient(rhoTable,obj.Tstepsize,obj.Pstepsize);
%             drhodp = interp2(Tgrid, Pgrid, table_drhodp, T, p, 'linear');
%         end
%         function [drhodT,d2rhodT2] = ComputeDrhoDT(obj, Pgrid, Tgrid, rhoTable, T, p)
%             [table_drhodT,~] = gradient(rhoTable,obj.Tstepsize,obj.Pstepsize);
%             drhodT = interp2(Tgrid, Pgrid, table_drhodT, T, p, 'linear');
%             [table_d2rhodT2,~] = gradient(table_drhodT,obj.Tstepsize,obj.Pstepsize); 
%             d2rhodT2 = interp2(Tgrid, Pgrid, table_d2rhodT2, T, p, 'linear');   
%         end
%         function dhdp = ComputeDhDp(obj, Pgrid, Tgrid, PhaseEnthalpyTable, T, p)
%             [~,table_dhdp] = gradient(PhaseEnthalpyTable,obj.Tstepsize,obj.Pstepsize);
%             dhdp = interp2(Tgrid, Pgrid, table_dhdp, T, p);
%         end
%         function [dhdT,d2hdT2] = ComputeDhDT(obj, Pgrid, Tgrid, PhaseEnthalpyTable, T, p)
%             [table_dhdT,~] = gradient(PhaseEnthalpyTable,obj.Tstepsize,obj.Pstepsize);
%             dhdT = interp2(Tgrid, Pgrid, table_dhdT, T, p, 'linear');
%             [table_d2hdT2,~] = gradient(table_dhdT,obj.Tstepsize,obj.Pstepsize); 
%             d2hdT2 = interp2(Tgrid, Pgrid, table_d2hdT2, T, p, 'linear');   
%         end
        
        
        % Injection properties; we are injecting only water, so it is
        % easier to use existing functions from Geothermal singlephase
        function rho = ComputeWaterDensity(obj, p, T)
            cw =  (0.0839*T.^2 - 64.593*T + 12437)*1e-12;
            rhofs = -0.0032*T.^2 + 1.7508*T + 757.5; 
            rho = rhofs.*(1+cw.*(p-obj.Psat));
        end
        function h = ComputeWaterEnthalpy(obj, p, T)
            rho = obj.ComputeWaterDensity(p, T);
            h = obj.uws + obj.Cp*(T-obj.Tsat)+p./rho;
        end
        function mu = ComputeWaterViscosity(obj, T)
            A = 2.414e-5;   B = 247.8;  C = T-140;   D = B./C;   E = 10.^D;            
            mu = A.*E;
        end

%         function rho = GetDensity(obj, Pgrid, Hgrid, rhoTable, h, p)
%             rho = interp2(Hgrid, Pgrid, rhoTable, h, p, 'linear');
% %             rho = interp2(Hgrid, Pgrid, rhoTable, h, p, 'linear');
% %             OTHER OPTIONS: 'makima' (gives NaN values...), 'cubic'
%         end
%         function mu = GetViscosity(obj, Pgrid, Hgrid, muTable, h, p)
%             mu = interp2(Hgrid, Pgrid, muTable, h, p, 'linear');
%         end

        function v = ComputeVelocity(obj, p, mu)
            % virtual call
        end
    end
end