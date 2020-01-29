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
        Kf % Conductivity
    end
    % the first derivative is an AVERAGE derivative
    methods
        function [S_w, S_s] = ComputeSaturation(obj, p, h) 
            S_w = TablePH.WaterSaturation(sub2ind(size(TablePH.WaterSaturation), obj.Pindex, obj.Hindex));
            S_s = TablePH.SteamSaturation(sub2ind(size(TablePH.SteamSaturation), obj.Pindex, obj.Hindex));        
        end
        % do we need dSdp_w and dSdh_w etc??
        
        function cond = AddConductivity(obj, p, T)
            cond = obj.Kf * ones(size(p));
        end
        
        function T = ComputeTemperature(obj, p, h)
            T = TablePH.Temperature(sub2ind(size(TablePH.Temperature), obj.Pindex, obj.Hindex));     
        end
        function [dTdp, d2Td2p] = ComputeDTDp(obj, p, h)
            % 1st derivative 
            [~,table_dTdp] = gradient(TablePH.Temperature,1,0.1); 
            dTdp = table_dTdp(sub2ind(size(table_dTdp), obj.Pindex, obj.Hindex));
            
            % 2nd derivative
            [~,table_d2Td2p] = gradient(table_dTdp,0.1); 
            d2Td2p = table_d2Td2p(sub2ind(size(table_d2Td2p), obj.Pindex, obj.Hindex));   
        end
        function [dTdh, d2Td2h] = ComputeDTDh(obj, p, h)
            % 1st derivative
            [table_dTdh,~] = gradient(TablePH.WaterDensity,1,0.1); %gradient('matrix','stepsize hor(j)','stepsize vert(i)') and you have P,H as i,j
            dTdh = table_dTdh(sub2ind(size(table_dTdh), obj.Pindex, obj.Hindex));
            
            % 2nd derivative
            [table_d2Td2h,~] = gradient(table_dTdh,1); %gradient('matrix','stepsize hor(j)','stepsize vert(i)') and you have P,H as i,j
            d2Td2h = table_d2Td2h(sub2ind(size(table_d2Td2h), obj.Pindex, obj.Hindex)); 
        end
        
        
        
        function rho = ComputeDensity(obj, Pindex, Hindex, rhoTable) % CORRECT !!
            rho = rhoTable( sub2ind(size(rhoTable), Pindex, Hindex) );
        end
        function [drhodp, d2rhod2p] = ComputeDrhoDp(obj, Pindex, Hindex, rhoTable) % CORRECT !!
            % derivative using matlab gradient() -> you must specify the step size!!
            [~,table_drhodp] = gradient(rhoTable,1,0.1); %gradient('matrix','stepsize hor(j)','stepsize vert(i)') and you have P,H as i,j
            drhodp = table_drhodp(sub2ind(size(table_drhodp), Pindex, Hindex));
            
            % 2nd derivatives
            [~,table_d2rhod2p] = gradient(table_drhodp,0.1); % specify stepsize for pressure (make this generic)
            d2rhod2p = table_d2rhod2p(sub2ind(size(table_d2rhod2p), Pindex, Hindex));
        end
        function [drhodh, d2rhod2h] = ComputeDrhoDh(obj, Pindex, Hindex, rhoTable) % CORRECT !!
            % 1st derivatives
            [table_drhodh,~] = gradient(rhoTable,1,0.1); %gradient('matrix','stepsize hor(j)','stepsize vert(i)') and you have P,H as i,j
            drhodh = table_drhodh(sub2ind(size(table_drhodh), Pindex, Hindex));
            
            % 2nd derivatives
            [table_d2rhod2h,~] = gradient(table_drhodh,1); % specify stepsize for enthalpy (make this generic)
            d2rhod2h = table_d2rhod2h(sub2ind(size(table_d2rhod2h), Pindex, Hindex));
        end
        
        
        
        
        
        function [mu_w, mu_s] = ComputeViscosity(obj, p, h)
            mu_w = TablePH.WaterViscosity(sub2ind(size(TablePH.WaterViscosity), obj.Pindex, obj.Hindex)); 
            mu_s = TablePH.SteamViscosity(sub2ind(size(TablePH.SteamViscosity), obj.Pindex, obj.Hindex)); 
        end
        function [dmudp_w, dmudp_s, d2mud2p_w, d2mud2p_s] = ComputeDmuDp(obj, p, h)
            % derivative using matlab gradient() -> you must specify the step size!!
            [~,table_dmudp_w] = gradient(TablePH.WaterViscosity,1,0.1); 
            [~,table_dmudp_s] = gradient(TablePH.SteamViscosity,1,0.1); 
            
            dmudp_w = table_dmudp_w(sub2ind(size(table_dmudp_w), obj.Pindex, obj.Hindex));
            dmudp_s = table_dmudp_s(sub2ind(size(table_dmudp_s), obj.Pindex, obj.Hindex));
            
            % 2nd derivatives
            [~,table_d2mud2p_w] = gradient(table_dmudp_w,0.1); 
            [~,table_d2mud2p_s] = gradient(table_dmudp_s,0.1);

            d2mud2p_w = table_d2mud2p_w(sub2ind(size(table_d2mud2p_w), obj.Pindex, obj.Hindex));
            d2mud2p_s = table_d2mud2p_s(sub2ind(size(table_d2mud2p_s), obj.Pindex, obj.Hindex));    
        end            
        function [dmudh_w, dmudh_s, d2mud2h_w, d2mud2h_s] = ComputeDmuDh(obj, p, h) 
            % 1st derivatives
            [table_dmudh_w,~] = gradient(TablePH.WaterViscosity,1,0.1); 
            [table_dmudh_s,~] = gradient(TablePH.SteamViscosity,1,0.1); 

            dmudh_w = table_dmudh_w(sub2ind(size(table_dmudh_w), obj.Pindex, obj.Hindex));
            dmudh_s = table_dmudh_s(sub2ind(size(table_dmudh_s), obj.Pindex, obj.Hindex));
            
            % 2nd derivatives
            [table_d2mud2h_w,~] = gradient(table_dmudh_w,1); 
            [table_d2mud2h_s,~] = gradient(table_dmudh_s,1);

            d2mud2h_w = table_d2mud2h_w(sub2ind(size(table_d2mud2h_w), obj.Pindex, obj.Hindex));
            d2mud2h_s = table_d2mud2h_s(sub2ind(size(table_d2mud2h_s), obj.Pindex, obj.Hindex));
        end
        
        function [U_w, U_s] = ComputeInternalEnergy(obj, p, h)
            U_w = TablePH.WaterInternalEnergy(sub2ind(size(TablePH.WaterInternalEnergy), obj.Pindex, obj.Hindex)); 
            U_s = TablePH.SteamInternalEnergy(sub2ind(size(TablePH.SteamInternalEnergy), obj.Pindex, obj.Hindex)); 
        end   
        function [dUdp_w, dUdp_s, d2Ud2p_w, d2Ud2p_s] = ComputeDUDp(obj, p, h)
            % 1st derivative 
            [~,table_dUdp_w] = gradient(TablePH.WaterInternalEnergy,1,0.1); 
            [~,table_dUdp_s] = gradient(TablePH.SteamInternalEnergy,1,0.1); 

            dUdp_w = table_dUdp_w(sub2ind(size(table_dUdp_w), obj.Pindex, obj.Hindex));
            dUdp_s = table_dUdp_s(sub2ind(size(table_dUdp_s), obj.Pindex, obj.Hindex));
            
            % 2nd derivatives
            [~,table_d2Ud2p_w] = gradient(table_dUdp_w,0.1); 
            [~,table_d2Ud2p_s] = gradient(table_dUdp_s,0.1);

            d2Ud2p_w = table_d2Ud2p_w(sub2ind(size(table_d2Ud2p_w), obj.Pindex, obj.Hindex));
            d2Ud2p_s = table_d2Ud2p_s(sub2ind(size(table_d2Ud2p_s), obj.Pindex, obj.Hindex));    
        end
        function [dUdh_w, dUdh_s, d2Ud2h_w, d2Ud2h_s] = ComputeDUDh(obj, p, h)
            % 1st derivatives
            [table_dUdh_w,~] = gradient(TablePH.WaterInternalEnergy,1,0.1); 
            [table_dUdh_s,~] = gradient(TablePH.SteamInternalEnergy,1,0.1); 

            dUdh_w = table_dUdh_w(sub2ind(size(table_dUdh_w), obj.Pindex, obj.Hindex));
            dUdh_s = table_dUdh_s(sub2ind(size(table_dUdh_s), obj.Pindex, obj.Hindex));
            
            % 2nd derivatives
            [table_d2Ud2h_w,~] = gradient(table_dUdh_w,1); 
            [table_d2Ud2h_s,~] = gradient(table_dUdh_s,1);

            d2Ud2h_w = table_d2Ud2h_w(sub2ind(size(table_d2Ud2h_w), obj.Pindex, obj.Hindex));
            d2Ud2h_s = table_d2Ud2h_s(sub2ind(size(table_d2Ud2h_s), obj.Pindex, obj.Hindex));
        end 
        
        function v = ComputeVelocity(obj, p, mu)
            % virtual call
        end
    end
end
