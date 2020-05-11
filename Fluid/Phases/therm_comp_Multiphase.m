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
%         Pstepsize = 1e4;
%         Hstepsize = 1e3; % implement this
        Pstepsize = 1e5;
        Hstepsize = 1e4; % implement this
        
        % for the injection well properties
        Cp_std              % Specific Heat of Phase in standard condition
        uws = 420000;         % internal energy at saturation J/kg
        Tsat = 373;         % T at saturation condition (assumed constant 100 C)
        Psat = 1e5;         % P at saturation condition (assumed 1e5 Pa)
    end
    % the first derivative is an AVERAGE derivative
    methods
        % Get Phase properties from tables ( LOOK_UP )
% %         function rho = GetDensity(obj, Pindex, Hindex, rhoTable)
% %             rho = rhoTable( sub2ind(size(rhoTable), Pindex, Hindex) );
% %         end
% %         function S = GetSaturation(obj, Pindex, Hindex, STable)
% %             S = STable( sub2ind(size(STable), Pindex, Hindex) );
% %         end 
% %         function mu = GetViscosity(obj, Pindex, Hindex, muTable)
% %             mu = muTable( sub2ind(size(muTable), Pindex, Hindex) );
% %         end
% %         function ThermCond = GetConductivity(obj, Pindex, Hindex, ThermCondTable)
% %             ThermCond = ThermCondTable( sub2ind(size(ThermCondTable), Pindex, Hindex));
% %         end
% %         function PhaseEnthalpy = GetPhaseEnthalpy(obj, Pindex, PhaseEnthalpyTable)
% %             PhaseEnthalpy = PhaseEnthalpyTable( sub2ind(size(PhaseEnthalpyTable), Pindex))'; %transpose to get column vector; is due to sub2ind thingy
% %         end
% %         
% %         % Derivatives (directly from tables)
% %         function drhodp = ComputeDrhoDp(obj, Pindex, Hindex, rhoTable)
% %             [~,table_drhodp] = gradient(rhoTable,obj.Hstepsize,obj.Pstepsize); 
% %             drhodp = table_drhodp(sub2ind(size(table_drhodp), Pindex, Hindex));      
% %         end
% %         function drhodh = ComputeDrhoDh(obj, Pindex, Hindex, rhoTable)
% %             [table_drhodh,~] = gradient(rhoTable,obj.Hstepsize,obj.Pstepsize); 
% %             drhodh = table_drhodh(sub2ind(size(table_drhodh), Pindex, Hindex));
% %         end
% %         function drho_times_Sdp = ComputeDrho_times_SDp(obj, Pindex, Hindex, rho_times_STable)
% %             [~,table_drho_times_Sdp] = gradient(rho_times_STable,obj.Hstepsize,obj.Pstepsize); 
% %             drho_times_Sdp = table_drho_times_Sdp(sub2ind(size(table_drho_times_Sdp), Pindex, Hindex));      
% %         end
% %         function drho_times_Sdh = ComputeDrho_times_SDh(obj, Pindex, Hindex, rho_times_STable) 
% %             [table_drho_times_Sdh,~] = gradient(rho_times_STable,obj.Hstepsize,obj.Pstepsize); 
% %             drho_times_Sdh = table_drho_times_Sdh(sub2ind(size(table_drho_times_Sdh), Pindex, Hindex));
% %         end
% %         function drho_times_hdp = ComputeDrho_times_hDp(obj, Pindex, Hindex, rho_times_hTable)
% %             [~,table_drho_times_hdp] = gradient(rho_times_hTable,obj.Hstepsize,obj.Pstepsize);
% %             drho_times_hdp = table_drho_times_hdp(sub2ind(size(table_drho_times_hdp), Pindex, Hindex));
% %         end
% %         function drho_times_hdh = ComputeDrho_times_hDh(obj, Pindex, Hindex, rho_times_hTable)
% %             [table_drho_times_hdh,~] = gradient(rho_times_hTable,obj.Hstepsize,obj.Pstepsize);
% %             drho_times_hdh = table_drho_times_hdh(sub2ind(size(table_drho_times_hdh), Pindex, Hindex));
% %         end
% %         function drhoHSdp = ComputeDrhoHSDp(obj, Pindex, Hindex, rhoHSTable)
% %             [~,table_drhoHSdp] = gradient(rhoHSTable,obj.Hstepsize,obj.Pstepsize);
% %             drhoHSdp = table_drhoHSdp(sub2ind(size(table_drhoHSdp), Pindex, Hindex));
% %         end
% %         function drhoHSdh = ComputeDrhoHSDh(obj, Pindex, Hindex, rhoHSTable)
% %             [table_drhoHSdh,~] = gradient(rhoHSTable,obj.Hstepsize,obj.Pstepsize);
% %             drhoHSdh = table_drhoHSdh(sub2ind(size(table_drhoHSdh), Pindex, Hindex));
% %         end
% %         function dTdh = ComputeDTDh(obj, Pindex, Hindex, TTable)
% %             [table_dTdh,~] = gradient(TTable,obj.Hstepsize,obj.Pstepsize); 
% %             dTdh = table_dTdh(sub2ind(size(table_dTdh), Pindex, Hindex));
% %         end
% %         function dTdp = ComputeDTDp(obj, Pindex, Hindex, TTable) 
% %             [~,table_dTdp] = gradient(TTable,obj.Hstepsize,obj.Pstepsize); 
% %             dTdp = table_dTdp(sub2ind(size(table_dTdp), Pindex, Hindex));
% %         end
% %         
% %         % i have copy pasted the newest derivatives from pressure to
% %         % enthalpy, so therefore the derivative with respect to enthalpy
% %         % was in the wrong position....
% %         function dmudp = ComputeDmuDp(obj, Pindex, Hindex, muTable)
% %             [~,table_dmudp] = gradient(muTable,obj.Hstepsize,obj.Pstepsize);
% %             dmudp = table_dmudp(sub2ind(size(table_dmudp), Pindex, Hindex));
% %         end       
% %         function dmudh = ComputeDmuDh(obj, Pindex, Hindex, muTable)
% %             [table_dmudh,~] = gradient(muTable,obj.Hstepsize,obj.Pstepsize);
% %             dmudh = table_dmudh(sub2ind(size(table_dmudh), Pindex, Hindex));
% %         end
% %         function dSdp = ComputeDSDp(obj, Pindex, Hindex, STable)
% %             [~,table_dSdp] = gradient(STable,obj.Hstepsize,obj.Pstepsize);
% %             dSdp = table_dSdp(sub2ind(size(table_dSdp), Pindex, Hindex));
% %         end
% %         function dSdh = ComputeDSDh(obj, Pindex, Hindex, STable)
% %             [table_dSdh,~] = gradient(STable,obj.Hstepsize,obj.Pstepsize);
% %             dSdh = table_dSdh(sub2ind(size(table_dSdh), Pindex, Hindex));
% %         end
% % 
% %         % These depend on how we treat the conductive flux term
% %         function d2Td2p = ComputeD2TD2p(obj, Pindex, Hindex, TTable) 
% %             [~,table_dTdp] = gradient(TTable,obj.Hstepsize,obj.Pstepsize); 
% %             [~,table_d2Td2p] = gradient(table_dTdp,obj.Hstepsize,obj.Pstepsize); 
% %             d2Td2p = table_d2Td2p(sub2ind(size(table_d2Td2p), Pindex, Hindex));   
% %         end
% %         function d2Td2h = ComputeD2TD2h(obj, Pindex, Hindex, TTable)
% %             [table_dTdh,~] = gradient(TTable,obj.Hstepsize,obj.Pstepsize); 
% %             [table_d2Td2h,~] = gradient(table_dTdh,obj.Hstepsize,obj.Pstepsize); 
% %             d2Td2h = table_d2Td2h(sub2ind(size(table_d2Td2h), Pindex, Hindex)); 
% %         end
% %         
% %         % 2nd derivatives for inflexion point correction
% %         function d2rhodp2 = ComputeD2rhoDp2(obj, Pindex, Hindex, rhoTable)
% %             % 1st derivative
% %             [~,table_drhodp] = gradient(rhoTable,obj.Hstepsize,obj.Pstepsize); 
% %             % 2nd derivative
% %             [~,table_d2rhodp2] = gradient(table_drhodp,obj.Hstepsize,obj.Pstepsize); % specify stepsize for pressure (make this generic)
% %             d2rhodp2 = table_d2rhodp2(sub2ind(size(table_d2rhodp2), Pindex, Hindex));
% %         end
% %         function d2rhodh2 = ComputeD2rhoDh2(obj, Pindex, Hindex, rhoTable)
% %             % 1st derivative
% %             [table_drhodh,~] = gradient(rhoTable,obj.Hstepsize,obj.Pstepsize); 
% %             % 2nd derivative
% %             [table_d2rhodh2,~] = gradient(table_drhodh,obj.Hstepsize,obj.Pstepsize); % specify stepsize for enthalpy (make this generic)
% %             d2rhodh2 = table_d2rhodh2(sub2ind(size(table_d2rhodh2), Pindex, Hindex));
% %         end
% %         function d2mudh2 = ComputeD2muDh2(obj, Pindex, Hindex, muTable)
% %             [table_dmudh,~] = gradient(muTable,obj.Hstepsize,obj.Pstepsize); 
% %             [table_d2mudh2,~] = gradient(table_dmudh,1); 
% %             d2mudh2 = table_d2mudh2(sub2ind(size(table_d2mudh2), Pindex, Hindex));
% %         end

        

%         Get properties using LINEAR INTERPOLATION
        function rho = GetDensity(obj, Pgrid, Hgrid, rhoTable, h, p)
            rho = interp2(Hgrid, Pgrid, rhoTable, h, p);
        end
        function S = GetSaturation(obj, Pgrid, Hgrid, STable, h, p)
            S = interp2(Hgrid, Pgrid, STable, h, p);
        end 
        function mu = GetViscosity(obj, Pgrid, Hgrid, muTable, h, p)
            mu = interp2(Hgrid, Pgrid, muTable, h, p);
        end
        function ThermCond = GetConductivity(obj, Pgrid, Hgrid, ThermCondTable, h, p)
            ThermCond = interp2(Hgrid, Pgrid, ThermCondTable, h, p);
        end
        function PhaseEnthalpy = GetPhaseEnthalpy(obj, Ptable, PhaseEnthalpyTable, p)
            PhaseEnthalpy = interp1(Ptable, PhaseEnthalpyTable' ,p); %transpose to get column vector; is due to sub2ind thingy
        end

        % Derivatives (directly from tables)
        function drhodp = ComputeDrhoDp(obj, Pgrid, Hgrid, rhoTable, h, p)
            [~,table_drhodp] = gradient(rhoTable,obj.Hstepsize,obj.Pstepsize); 
            drhodp = interp2(Hgrid, Pgrid, table_drhodp, h, p);
        end
        function drhodh = ComputeDrhoDh(obj, Pgrid, Hgrid, rhoTable, h, p)
            [table_drhodh,~] = gradient(rhoTable,obj.Hstepsize,obj.Pstepsize); 
            drhodh = interp2(Hgrid, Pgrid, table_drhodh, h, p);
        end
        function drho_times_Sdp = ComputeDrho_times_SDp(obj, Pgrid, Hgrid, rho_times_STable, h, p)
            [~,table_drho_times_Sdp] = gradient(rho_times_STable,obj.Hstepsize,obj.Pstepsize); 
            drho_times_Sdp = interp2(Hgrid, Pgrid, table_drho_times_Sdp, h, p);     
        end
        function drho_times_Sdh = ComputeDrho_times_SDh(obj, Pgrid, Hgrid, rho_times_STable, h, p) 
            [table_drho_times_Sdh,~] = gradient(rho_times_STable,obj.Hstepsize,obj.Pstepsize); 
            drho_times_Sdh = interp2(Hgrid, Pgrid, table_drho_times_Sdh, h, p);
        end
        function drho_times_hdp = ComputeDrho_times_hDp(obj, Pgrid, Hgrid, rho_times_hTable, h, p)
            [~,table_drho_times_hdp] = gradient(rho_times_hTable,obj.Hstepsize,obj.Pstepsize);
            drho_times_hdp = interp2(Hgrid, Pgrid, table_drho_times_hdp, h, p);
        end
        function drho_times_hdh = ComputeDrho_times_hDh(obj, Pgrid, Hgrid, rho_times_hTable, h, p)
            [table_drho_times_hdh,~] = gradient(rho_times_hTable,obj.Hstepsize,obj.Pstepsize);
            drho_times_hdh = interp2(Hgrid, Pgrid, table_drho_times_hdh, h, p);
        end
        function drhoHSdp = ComputeDrhoHSDp(obj, Pgrid, Hgrid, rhoHSTable, h, p)
            [~,table_drhoHSdp] = gradient(rhoHSTable,obj.Hstepsize,obj.Pstepsize);
            drhoHSdp = interp2(Hgrid, Pgrid, table_drhoHSdp, h, p);
        end
        function drhoHSdh = ComputeDrhoHSDh(obj, Pgrid, Hgrid, rhoHSTable, h, p)
            [table_drhoHSdh,~] = gradient(rhoHSTable,obj.Hstepsize,obj.Pstepsize);
            drhoHSdh = interp2(Hgrid, Pgrid, table_drhoHSdh, h, p);
        end
        function dTdh = ComputeDTDh(obj, Pgrid, Hgrid, TTable, h, p)
            [table_dTdh,~] = gradient(TTable,obj.Hstepsize,obj.Pstepsize); 
            dTdh = interp2(Hgrid, Pgrid, table_dTdh, h, p);
        end
        function dTdp = ComputeDTDp(obj, Pgrid, Hgrid, TTable, h, p) 
            [~,table_dTdp] = gradient(TTable,obj.Hstepsize,obj.Pstepsize); 
            dTdp = interp2(Hgrid, Pgrid, table_dTdp, h, p);
        end
        
        % i have copy pasted the newest derivatives from pressure to
        % enthalpy, so therefore the derivative with respect to enthalpy
        % was in the wrong position....
        function dmudp = ComputeDmuDp(obj, Pgrid, Hgrid, muTable, h, p)
            [~,table_dmudp] = gradient(muTable,obj.Hstepsize,obj.Pstepsize);
            dmudp = interp2(Hgrid, Pgrid, table_dmudp, h, p);
        end       
        function dmudh = ComputeDmuDh(obj, Pgrid, Hgrid, muTable, h, p)
            [table_dmudh,~] = gradient(muTable,obj.Hstepsize,obj.Pstepsize);
            dmudh = interp2(Hgrid, Pgrid, table_dmudh, h, p);
        end
        function dSdp = ComputeDSDp(obj, Pgrid, Hgrid, STable, h, p)
            [~,table_dSdp] = gradient(STable,obj.Hstepsize,obj.Pstepsize);
            dSdp = interp2(Hgrid, Pgrid, table_dSdp, h, p);
        end
        function dSdh = ComputeDSDh(obj, Pgrid, Hgrid, STable, h, p)
            [table_dSdh,~] = gradient(STable,obj.Hstepsize,obj.Pstepsize);
            dSdh = interp2(Hgrid, Pgrid, table_dSdh, h, p);
        end

        % These depend on how we treat the conductive flux term
        function d2Td2p = ComputeD2TD2p(obj, Pgrid, Hgrid, TTable, h, p) 
            [~,table_dTdp] = gradient(TTable,obj.Hstepsize,obj.Pstepsize); 
            [~,table_d2Td2p] = gradient(table_dTdp,obj.Hstepsize,obj.Pstepsize); 
            d2Td2p = interp2(Hgrid, Pgrid, table_d2Td2p, h, p);   
        end
        function d2Td2h = ComputeD2TD2h(obj, Pgrid, Hgrid, TTable, h, p)
            [table_dTdh,~] = gradient(TTable,obj.Hstepsize,obj.Pstepsize); 
            [table_d2Td2h,~] = gradient(table_dTdh,obj.Hstepsize,obj.Pstepsize); 
            d2Td2h = interp2(Hgrid, Pgrid, table_d2Td2h, h, p); 
        end
        
        % 2nd derivatives for inflexion point correction
        function d2rhodp2 = ComputeD2rhoDp2(obj, Pgrid, Hgrid, rhoTable, h, p)
            % 1st derivative
            [~,table_drhodp] = gradient(rhoTable,obj.Hstepsize,obj.Pstepsize); 
            % 2nd derivative
            [~,table_d2rhodp2] = gradient(table_drhodp,obj.Hstepsize,obj.Pstepsize); % specify stepsize for pressure (make this generic)
            d2rhodp2 = interp2(Hgrid, Pgrid, table_d2rhodp2, h, p);
        end
        function d2rhodh2 = ComputeD2rhoDh2(obj, Pgrid, Hgrid, rhoTable, h, p)
            % 1st derivative
            [table_drhodh,~] = gradient(rhoTable,obj.Hstepsize,obj.Pstepsize); 
            % 2nd derivative
            [table_d2rhodh2,~] = gradient(table_drhodh,obj.Hstepsize,obj.Pstepsize); % specify stepsize for enthalpy (make this generic)
            d2rhodh2 = interp2(Hgrid, Pgrid, table_d2rhodh2, h, p);
        end
        function d2mudh2 = ComputeD2muDh2(obj, Pgrid, Hgrid, muTable, h, p)
            [table_dmudh,~] = gradient(muTable,obj.Hstepsize,obj.Pstepsize); 
            [table_d2mudh2,~] = gradient(table_dmudh,1); 
            d2mudh2 = interp2(Hgrid, Pgrid, table_d2mudh2, h, p);
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
%         function obj = ComputeEnthalpy(obj, p, T)
%             % virtual call
%         end

        
        
        
        
        %%% ORIGINALLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function U = GetInternalEnergy(obj, Pindex, Hindex, UTable)
            U = UTable( sub2ind(size(UTable), Pindex, Hindex) );
        end
        function drho_over_mudp = ComputeDrho_over_muDp(obj, Pindex, Hindex, rho_over_muTable)
            [~,table_drho_over_mudp] = gradient(rho_over_muTable,obj.Pstepsize);
            drho_over_mudp = table_drho_over_mudp(sub2ind(size(table_drho_over_mudp), Pindex, Hindex));
        end       
        function dhdp = ComputeDhDp(obj, Pindex, Hindex, hTable)
            [~,table_dhdp] = gradient(hTable,obj.Pstepsize);
            dhdp = table_dhdp(sub2ind(size(table_dhdp), Pindex, Hindex));
        end
        function dUfdp = ComputeDUfDp(obj, Pindex, Hindex, UfTable)
            [~,table_dUfdp] = gradient(UfTable,obj.Pstepsize);
            dUfdp = table_dUfdp(sub2ind(size(table_dUfdp), Pindex, Hindex));
        end
        function dUfdh = ComputeDUfDh(obj, Pindex, Hindex, UfTable)
            [table_dUfdh,~] = gradient(UfTable,obj.Hstepsize);
            dUfdh = table_dUfdh(sub2ind(size(table_dUfdh), Pindex, Hindex));
        end        
        % For the thermal conductivity tensor:
        function ds_times_conddp = ComputeDs_times_condDp(obj, Pindex, Hindex, s_times_condTable)
            [~,table_ds_times_conddp] = gradient(s_times_condTable,obj.Hstepsize,obj.Pstepsize);
            ds_times_conddp = table_ds_times_conddp(sub2ind(size(table_ds_times_conddp), Pindex, Hindex));
        end
        function ds_times_conddh = ComputeDs_times_condDh(obj, Pindex, Hindex, s_times_condTable)
            [table_ds_times_conddh,~] = gradient(s_times_condTable,obj.Hstepsize,obj.Pstepsize);
            ds_times_conddh = table_ds_times_conddh(sub2ind(size(table_ds_times_conddh), Pindex, Hindex));
        end
        
                        
                
    end
end
