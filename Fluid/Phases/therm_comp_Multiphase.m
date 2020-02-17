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
        function rho = ComputeDensity(obj, Pindex, Hindex, rhoTable)
            rho = rhoTable( sub2ind(size(rhoTable), Pindex, Hindex) );
        end
        function S = ComputeSaturation(obj, Pindex, Hindex, STable)
            S = STable( sub2ind(size(STable), Pindex, Hindex) );
        end 
        function mu = ComputeViscosity(obj, Pindex, Hindex, muTable)
            mu = muTable( sub2ind(size(muTable), Pindex, Hindex) );
        end
        function U = ComputeInternalEnergy(obj, Pindex, Hindex, UTable)
            U = UTable( sub2ind(size(UTable), Pindex, Hindex) );
        end   

        function cond = AddConductivity(obj, p, T)
            cond = obj.Kf * ones(size(p));
        end
        
        function [drhodp, d2rhod2p] = ComputeDrhoDp(obj, Pindex, Hindex, rhoTable)
            % derivative using matlab gradient() -> you must specify the step size!!
            [~,table_drhodp] = gradient(rhoTable,1,0.1); %gradient('matrix','stepsize hor(j)','stepsize vert(i)') and you have P,H as i,j
            drhodp = table_drhodp(sub2ind(size(table_drhodp), Pindex, Hindex));      
            % 2nd derivative
            [~,table_d2rhod2p] = gradient(table_drhodp,0.1); % specify stepsize for pressure (make this generic)
            d2rhod2p = table_d2rhod2p(sub2ind(size(table_d2rhod2p), Pindex, Hindex));
        end
        function [drhodh, d2rhod2h] = ComputeDrhoDh(obj, Pindex, Hindex, rhoTable)
            % 1st derivative
            [table_drhodh,~] = gradient(rhoTable,1,0.1); 
            drhodh = table_drhodh(sub2ind(size(table_drhodh), Pindex, Hindex));
            
            % 2nd derivative
            [table_d2rhod2h,~] = gradient(table_drhodh,1); % specify stepsize for enthalpy (make this generic)
            d2rhod2h = table_d2rhod2h(sub2ind(size(table_d2rhod2h), Pindex, Hindex));
        end
        function [dSdp, d2Sd2p] = ComputeDSDp(obj, Pindex, Hindex, STable)
            [~,table_dSdp] = gradient(STable,1,0.1); 
            dSdp = table_dSdp(sub2ind(size(table_dSdp), Pindex, Hindex));
            [~,table_d2Sd2p] = gradient(table_dSdp,0.1); 
            d2Sd2p = table_d2Sd2p(sub2ind(size(table_d2Sd2p), Pindex, Hindex));
        end
        function [dSdh, d2Sd2h] = ComputeDSDh(obj, Pindex, Hindex, STable)
            [~,table_dSdh] = gradient(STable,1,0.1); 
            dSdh = table_dSdh(sub2ind(size(table_dSdh), Pindex, Hindex));
            [~,table_d2Sd2h] = gradient(table_dSdh,1); 
            d2Sd2h = table_d2Sd2h(sub2ind(size(table_d2Sd2h), Pindex, Hindex));
        end
        function [dmudp, d2mud2p] = ComputeDmuDp(obj, Pindex, Hindex, muTable)
            [~,table_dmudp] = gradient(muTable,1,0.1);             
            dmudp = table_dmudp(sub2ind(size(table_dmudp), Pindex, Hindex));
            [~,table_d2mud2p] = gradient(table_dmudp,0.1); 
            d2mud2p = table_d2mud2p(sub2ind(size(table_d2mud2p), Pindex, Hindex));
        end            
        function [dmudh, d2mud2h] = ComputeDmuDh(obj, Pindex, Hindex, muTable)
            [table_dmudh,~] = gradient(muTable,1,0.1); 
            dmudh = table_dmudh(sub2ind(size(table_dmudh), Pindex, Hindex));
            [table_d2mud2h,~] = gradient(table_dmudh,1); 
            d2mud2h = table_d2mud2h(sub2ind(size(table_d2mud2h), Pindex, Hindex));
        end
        function [dUdp, d2Ud2p] = ComputeDUDp(obj, Pindex, Hindex, UTable)
            [~,table_dUdp] = gradient(UTable,1,0.1); 
            dUdp = table_dUdp(sub2ind(size(table_dUdp), Pindex, Hindex));
            [~,table_d2Ud2p] = gradient(table_dUdp,0.1); 
            d2Ud2p = table_d2Ud2p(sub2ind(size(table_d2Ud2p), Pindex, Hindex));
        end
        function [dUdh, d2Ud2h] = ComputeDUDh(obj, Pindex, Hindex, UTable)
            [table_dUdh,~] = gradient(UTable,1,0.1); 
            dUdh = table_dUdh(sub2ind(size(table_dUdh), Pindex, Hindex));
            [table_d2Ud2h,~] = gradient(table_dUdh,1); 
            d2Ud2h = table_d2Ud2h(sub2ind(size(table_d2Ud2h), Pindex, Hindex));
        end 
        % Compute Temperature derivatives
        function [dTdp, d2Td2p] = ComputeDTDp(obj, Pindex, Hindex, TTable) 
            [~,table_dTdp] = gradient(TTable,1,0.1); 
            dTdp = table_dTdp(sub2ind(size(table_dTdp), Pindex, Hindex));
            [~,table_d2Td2p] = gradient(table_dTdp,0.1); 
            d2Td2p = table_d2Td2p(sub2ind(size(table_d2Td2p), Pindex, Hindex));   
        end
        function [dTdh, d2Td2h] = ComputeDTDh(obj, Pindex, Hindex, TTable)
            [table_dTdh,~] = gradient(TTable,1,0.1); 
            dTdh = table_dTdh(sub2ind(size(table_dTdh), Pindex, Hindex));
            [table_d2Td2h,~] = gradient(table_dTdh,1); 
            d2Td2h = table_d2Td2h(sub2ind(size(table_d2Td2h), Pindex, Hindex)); 
        end
        
        function v = ComputeVelocity(obj, p, mu)
            % virtual call
        end
    end
end
