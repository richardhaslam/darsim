% Thermal compressible phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Arjan Marelis
%TU Delft
%Created: 21 January 2020
%Last modified: 23 January 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef therm_comp_Multiphase < phase
    % Here, you only place functions for phase properties
    properties
        Pstepsize = 1e5;
        Hstepsize = 1e4; % implement this
        
        % for the injection well properties
        Cp_std              % Specific Heat of Phase in standard condition
        uws = 420000;         % internal energy at saturation J/kg
        Tsat = 373;         % T at saturation condition (assumed constant 100 C)
        Psat = 1e5;         % P at saturation condition (assumed 1e5 Pa)
    end
    methods
        function rho = GetDensity(obj, Pgrid, Hgrid, rhoTable, h, p)
            rho = interp2(Hgrid, Pgrid, rhoTable, h, p, 'linear');
%             rho = interp2(Hgrid, Pgrid, rhoTable, h, p, 'linear');
%             OTHER OPTIONS: 'makima' (gives NaN values...), 'cubic'
        end
        function S = GetSaturation(obj, Pgrid, Hgrid, STable, h, p)
            S = interp2(Hgrid, Pgrid, STable, h, p, 'linear');
        end 
        function mu = GetViscosity(obj, Pgrid, Hgrid, muTable, h, p)
            mu = interp2(Hgrid, Pgrid, muTable, h, p, 'linear');
        end
        function ThermCond = GetConductivity(obj, Pgrid, Hgrid, ThermCondTable, h, p)
            ThermCond = interp2(Hgrid, Pgrid, ThermCondTable, h, p, 'linear');
        end
        function PhaseEnthalpy = GetPhaseEnthalpy(obj, Ptable, PhaseEnthalpyTable, p)
            PhaseEnthalpy = interp1(Ptable, PhaseEnthalpyTable' ,p); %transpose to get column vector; is due to sub2ind thingy
        end

        % Derivatives (directly from tables)
        function drhodp = ComputeDrhoDp(obj, Pgrid, Hgrid, rhoTable, h, p)
            [~,table_drhodp] = gradient(rhoTable,obj.Hstepsize,obj.Pstepsize); 
            drhodp = interp2(Hgrid, Pgrid, table_drhodp, h, p, 'linear');
        end
        function drhodh = ComputeDrhoDh(obj, Pgrid, Hgrid, rhoTable, h, p)
            [table_drhodh,~] = gradient(rhoTable,obj.Hstepsize,obj.Pstepsize); 
            drhodh = interp2(Hgrid, Pgrid, table_drhodh, h, p, 'linear');
        end
        function drho_times_Sdp = ComputeDrho_times_SDp(obj, Pgrid, Hgrid, rho_times_STable, h, p)
            [~,table_drho_times_Sdp] = gradient(rho_times_STable,obj.Hstepsize,obj.Pstepsize); 
            drho_times_Sdp = interp2(Hgrid, Pgrid, table_drho_times_Sdp, h, p, 'linear');     
        end
        function drho_times_Sdh = ComputeDrho_times_SDh(obj, Pgrid, Hgrid, rho_times_STable, h, p) 
            [table_drho_times_Sdh,~] = gradient(rho_times_STable,obj.Hstepsize,obj.Pstepsize); 
            drho_times_Sdh = interp2(Hgrid, Pgrid, table_drho_times_Sdh, h, p, 'linear');
        end
        function drho_times_hdp = ComputeDrho_times_hDp(obj, Pgrid, Hgrid, rho_times_hTable, h, p)
            [~,table_drho_times_hdp] = gradient(rho_times_hTable,obj.Hstepsize,obj.Pstepsize);
            drho_times_hdp = interp2(Hgrid, Pgrid, table_drho_times_hdp, h, p, 'linear');
        end
        function drho_times_hdh = ComputeDrho_times_hDh(obj, Pgrid, Hgrid, rho_times_hTable, h, p)
            [table_drho_times_hdh,~] = gradient(rho_times_hTable,obj.Hstepsize,obj.Pstepsize);
            drho_times_hdh = interp2(Hgrid, Pgrid, table_drho_times_hdh, h, p, 'linear');
        end
        function drhoHSdp = ComputeDrhoHSDp(obj, Pgrid, Hgrid, rhoHSTable, h, p)
            [~,table_drhoHSdp] = gradient(rhoHSTable,obj.Hstepsize,obj.Pstepsize);
            drhoHSdp = interp2(Hgrid, Pgrid, table_drhoHSdp, h, p, 'linear');
        end
        function drhoHSdh = ComputeDrhoHSDh(obj, Pgrid, Hgrid, rhoHSTable, h, p)
            [table_drhoHSdh,~] = gradient(rhoHSTable,obj.Hstepsize,obj.Pstepsize);
            drhoHSdh = interp2(Hgrid, Pgrid, table_drhoHSdh, h, p, 'linear');
        end
        function dTdh = ComputeDTDh(obj, Pgrid, Hgrid, TTable, h, p)
            [table_dTdh,~] = gradient(TTable,obj.Hstepsize,obj.Pstepsize); 
            dTdh = interp2(Hgrid, Pgrid, table_dTdh, h, p, 'linear');
        end
        function dTdp = ComputeDTDp(obj, Pgrid, Hgrid, TTable, h, p) 
            [~,table_dTdp] = gradient(TTable,obj.Hstepsize,obj.Pstepsize); 
            dTdp = interp2(Hgrid, Pgrid, table_dTdp, h, p, 'linear');
        end
        
        function dmudp = ComputeDmuDp(obj, Pgrid, Hgrid, muTable, h, p)
            [~,table_dmudp] = gradient(muTable,obj.Hstepsize,obj.Pstepsize);
            dmudp = interp2(Hgrid, Pgrid, table_dmudp, h, p, 'linear');
        end       
        function dmudh = ComputeDmuDh(obj, Pgrid, Hgrid, muTable, h, p)
            [table_dmudh,~] = gradient(muTable,obj.Hstepsize,obj.Pstepsize);
            dmudh = interp2(Hgrid, Pgrid, table_dmudh, h, p, 'linear');
        end
        function dSdp = ComputeDSDp(obj, Pgrid, Hgrid, STable, h, p)
            [~,table_dSdp] = gradient(STable,obj.Hstepsize,obj.Pstepsize);
            dSdp = interp2(Hgrid, Pgrid, table_dSdp, h, p, 'linear');
        end
        function dSdh = ComputeDSDh(obj, Pgrid, Hgrid, STable, h, p)
            [table_dSdh,~] = gradient(STable,obj.Hstepsize,obj.Pstepsize);
            dSdh = interp2(Hgrid, Pgrid, table_dSdh, h, p, 'linear');
        end

        % These depend on how we treat the conductive flux term
        function d2Td2p = ComputeD2TD2p(obj, Pgrid, Hgrid, TTable, h, p) 
            [~,table_dTdp] = gradient(TTable,obj.Hstepsize,obj.Pstepsize); 
            [~,table_d2Td2p] = gradient(table_dTdp,obj.Hstepsize,obj.Pstepsize); 
            d2Td2p = interp2(Hgrid, Pgrid, table_d2Td2p, h, p, 'linear');   
        end
        function d2Td2h = ComputeD2TD2h(obj, Pgrid, Hgrid, TTable, h, p)
            [table_dTdh,~] = gradient(TTable,obj.Hstepsize,obj.Pstepsize); 
            [table_d2Td2h,~] = gradient(table_dTdh,obj.Hstepsize,obj.Pstepsize); 
            d2Td2h = interp2(Hgrid, Pgrid, table_d2Td2h, h, p, 'linear'); 
        end
        
        % 2nd derivatives for inflexion point correction
        function d2rhodp2 = ComputeD2rhoDp2(obj, Pgrid, Hgrid, rhoTable, h, p)
            % 1st derivative
            [~,table_drhodp] = gradient(rhoTable,obj.Hstepsize,obj.Pstepsize); 
            % 2nd derivative
            [~,table_d2rhodp2] = gradient(table_drhodp,obj.Hstepsize,obj.Pstepsize); % specify stepsize for pressure (make this generic)
            d2rhodp2 = interp2(Hgrid, Pgrid, table_d2rhodp2, h, p, 'linear');
        end
        function d2rhodh2 = ComputeD2rhoDh2(obj, Pgrid, Hgrid, rhoTable, h, p)
            % 1st derivative
            [table_drhodh,~] = gradient(rhoTable,obj.Hstepsize,obj.Pstepsize); 
            % 2nd derivative
            [table_d2rhodh2,~] = gradient(table_drhodh,obj.Hstepsize,obj.Pstepsize); % specify stepsize for enthalpy (make this generic)
            d2rhodh2 = interp2(Hgrid, Pgrid, table_d2rhodh2, h, p, 'linear');
        end
        function d2mudh2 = ComputeD2muDh2(obj, Pgrid, Hgrid, muTable, h, p)
            [table_dmudh,~] = gradient(muTable,obj.Hstepsize,obj.Pstepsize); 
            [table_d2mudh2,~] = gradient(table_dmudh,1); 
            d2mudh2 = interp2(Hgrid, Pgrid, table_d2mudh2, h, p, 'linear');
        end
        
        function dhdp = ComputeDhDp(obj, Ptable, hTable, p)
            table_dhdp = gradient(hTable,obj.Pstepsize);
            dhdp = interp1(Ptable, table_dhdp, p);
            
%             dhdp = ((rho - p.*drhodp)./rho.^2);

        end

        
        % Injection properties; we are injecting only water, so it is
        % easier to use existing functions from Geothermal singlephase
        function rho = ComputeWaterDensity(obj, p, T)
            cw =  (0.0839*T.^2 - 64.593*T + 12437)*1e-12;
            rhofs = -0.0032*T.^2 + 1.7508*T + 757.5; 
            rho = rhofs.*(1+cw.*(p-obj.Psat));
        end
        function h = ComputeWaterEnthalpy(obj, p, T)
            rho = obj.ComputeWaterDensity(p, T);
            h = obj.uws + obj.Cp_std*(T-obj.Tsat)+p./rho;
        end
        function mu = ComputeWaterViscosity(obj, T)
            A = 2.414e-5;   B = 247.8;  C = T-140;   D = B./C;   E = 10.^D;            
            mu = A.*E;
        end

        
        % Other
        function v = ComputeVelocity(obj, p, mu)
            % virtual call
        end
        function obj = ComputeDensity(obj)
            % virtual call
        end
        
    end
end
