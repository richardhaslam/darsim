% Thermal compressible phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Arjan Marelis
%TU Delft
%Created: 21 January 2020
%Last modified: 23 January 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef therm_comp_Multiphase < phase
    properties
        Pstepsize = 1e4;
        Hstepsize = 1e3; % implement this
    end
    % the first derivative is an AVERAGE derivative
    methods
        function rho = GetDensity(obj, Pindex, Hindex, rhoTable)
            rho = rhoTable( sub2ind(size(rhoTable), Pindex, Hindex) );
        end
        function S = GetSaturation(obj, Pindex, Hindex, STable)
            S = STable( sub2ind(size(STable), Pindex, Hindex) );
        end 
        function mu = GetViscosity(obj, Pindex, Hindex, muTable)
            mu = muTable( sub2ind(size(muTable), Pindex, Hindex) );
        end
        function U = GetInternalEnergy(obj, Pindex, Hindex, UTable)
            U = UTable( sub2ind(size(UTable), Pindex, Hindex) );
        end   
        function ThermCond = GetConductivity(obj, Pindex, Hindex, ThermCondTable)
            ThermCond = ThermCondTable( sub2ind(size(ThermCondTable), Pindex, Hindex));
        end
        function PhaseEnthalpy = GetPhaseEnthalpy(obj, Pindex, Hindex, PhaseEnthalpyTable)
            PhaseEnthalpy = PhaseEnthalpyTable( sub2ind(size(PhaseEnthalpyTable), Pindex, Hindex));
        end
        
        % Derivatives
        function drhodp = ComputeDrhoDp(obj, Pindex, Hindex, rhoTable)
            [~,table_drhodp] = gradient(rhoTable,obj.Pstepsize); 
            drhodp = table_drhodp(sub2ind(size(table_drhodp), Pindex, Hindex));      
        end
        function drhodh = ComputeDrhoDh(obj, Pindex, Hindex, rhoTable)
            [table_drhodh,~] = gradient(rhoTable,obj.Hstepsize); 
            drhodh = table_drhodh(sub2ind(size(table_drhodh), Pindex, Hindex));
        end
        function drhoTdp = ComputeDrhoTDp(obj, Pindex, Hindex, rhoTTable)
            [~,table_drhoTdp] = gradient(rhoTTable,obj.Pstepsize); 
            drhoTdp = table_drhoTdp(sub2ind(size(table_drhoTdp), Pindex, Hindex));      
        end
        function drhoTdh = ComputeDrhoTDh(obj, Pindex, Hindex, rhoTTable)
            [table_drhoTdh,~] = gradient(rhoTTable,obj.Hstepsize); 
            drhoTdh = table_drhodh(sub2ind(size(table_drhoTdh), Pindex, Hindex));
        end
        function drho_over_mudp = ComputeDrho_over_muDp(obj, Pindex, Hindex, rho_over_muTable)
            [~,table_drho_over_mudp] = gradient(rho_over_muTable,obj.Pstepsize);
            drho_over_mudp = table_drho_over_mudp(sub2ind(size(table_drho_over_mudp), Pindex, Hindex));
        end        
        function drho_times_hdp = ComputeDrho_times_hDp(obj, Pindex, Hindex, rho_times_hTable)
            [~,table_drho_times_hdp] = gradient(rho_times_hTable,obj.Pstepsize);
            drho_times_hdp = table_drho_times_hdp(sub2ind(size(table_drho_times_hdp), Pindex, Hindex));
        end
        function drho_times_hdh = ComputeDrho_times_hDh(obj, Pindex, Hindex, rho_times_hTable)
            [~,table_drho_times_hdh] = gradient(rho_times_hTable,obj.Pstepsize);
            drho_times_hdh = table_drho_times_hdh(sub2ind(size(table_drho_times_hdh), Pindex, Hindex));
        end
        function dhdp = ComputeDhDp(obj, Pindex, Hindex, hTable)
            [~,table_dhdp] = gradient(hTable,obj.Pstepsize);
            dhdp = table_dhdp(sub2ind(size(table_dhdp), Pindex, Hindex));
        end        
        function dUfdp = ComputeDUfDp(obj, Pindex, Hindex, UfTable)
            [table_dUfdp,~] = gradient(UfTable,obj.Hstepsize); 
            dUfdp = table_dUfdp(sub2ind(size(table_dUfdp), Pindex, Hindex));
        end
        function dUfdh = ComputeDUfDh(obj, Pindex, Hindex, UfTable)
            [table_dUfdh,~] = gradient(UfTable,obj.Hstepsize); 
            dUfdh = table_dUfdh(sub2ind(size(table_dUfdh), Pindex, Hindex));
        end        
        function [dTdp, d2Td2p] = ComputeDTDp(obj, Pindex, Hindex, TTable) 
            [~,table_dTdp] = gradient(TTable,obj.Pstepsize); 
            dTdp = table_dTdp(sub2ind(size(table_dTdp), Pindex, Hindex));
            [~,table_d2Td2p] = gradient(table_dTdp,obj.Pstepsize); 
            d2Td2p = table_d2Td2p(sub2ind(size(table_d2Td2p), Pindex, Hindex));   
        end
        function [dTdh, d2Td2h] = ComputeDTDh(obj, Pindex, Hindex, TTable)
            [table_dTdh,~] = gradient(TTable,obj.Hstepsize); 
            dTdh = table_dTdh(sub2ind(size(table_dTdh), Pindex, Hindex));
            [table_d2Td2h,~] = gradient(table_dTdh,obj.Hstepsize); 
            d2Td2h = table_d2Td2h(sub2ind(size(table_d2Td2h), Pindex, Hindex)); 
        end
        % For the thermal conductivity tensor:
        function ds_times_conddp = ComputeDs_times_condDp(obj, Pindex, Hindex, s_times_condTable)
            [~,table_ds_times_conddp] = gradient(s_times_condTable,obj.Pstepsize);
            ds_times_conddp = table_ds_times_conddp(sub2ind(size(table_ds_times_conddp), Pindex, Hindex));
        end
        function ds_times_conddh = ComputeDs_times_condDh(obj, Pindex, Hindex, s_times_condTable)
            [~,table_ds_times_conddh] = gradient(s_times_condTable,obj.Pstepsize);
            ds_times_conddh = table_ds_times_conddh(sub2ind(size(table_ds_times_conddh), Pindex, Hindex));
        end

        
        function v = ComputeVelocity(obj, p, mu)
            % virtual call
        end
    end
end
