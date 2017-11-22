% Brooks-Corey relative permeability model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Jeroen Rijntjes
%TU Delft
%Created: 16 October 2017
%Last modified: 23 October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef relperm_model_brookscorey < relperm_model
    properties
        n_w = 4;
        n_nw = 2; 
        lambda = 3;
    end
        
    methods
        
        function kr = ComputeRelPerm(obj, Phases, s)     
            % Wetting phase relative permeability            
            kr(:,1) = ((s-Phases(2).sr)/(1-Phases(2).sr)).^(obj.n_w);           
            kr(s < Phases(2).sr, 1) = 0;
           
            % Non-Wetting phase relative permeability
            kr(:,2) = ((1-s)/(1-Phases(2).sr)).^2 .* ( 1 - ((s-Phases(2).sr)/(1-Phases(2).sr)).^(obj.n_nw) );
            %kr(s > 1 - Phases(1).sr, 2) = 0;
                       
            kr(s < Phases(2).sr, 2) = 1;
            kr(:,2) = smooth(kr(:,2));
            
            %kr(s > 1 - Phases(1).sr,1) = 1;
            
%             if s > 1 - Phases(1).sr;
%             IDX = find(s == 1 - Phases(1).sr);
%             X = 1001 - IDX;
%             ValueKr = find (Kr(IDX,1));
%             Kr_add=linspace(ValueKr, 1, X); 
%             end
            
            
%             if s == 1 - Phases(1).sr;
%                new = find(s(1,:)>1-Phases(1).sr);
 
%                new = last row kr column; getindex; 
%                new = new + 1;
%                kr(new,1) = kr (previous) + 0.01;
%                new = new +1;
%                finish when kr (...,...) = 1
%             end  
%             end    
                         
            % lambda formulations                  
%             n_x = ((2+3*obj.lambda) / obj.lambda);
%             n_y = ((2+obj.lambda) / obj.lambda);
%             kr(:,3) = ((s-Phases(1).sr)/(1-Phases(1).sr)).^(n_x);
%             kr(:,4) = ((1-s)/(1-Phases(1).sr)).^2 .* ( 1 - ((s-Phases(1).sr)/(1-Phases(1).sr)).^(n_y) );
%             kr(s < Phases(1).sr, 3) = 0;
%             kr(s > 1 - Phases(2).sr, 4) = 0;
%             kr(s < Phases(1).sr, 4) = 1;
%             kr(s > 1 - Phases(2).sr, 3) = 1;
            
            
        end
        function dkr = ComputeDerivative(obj, Phases, s)
%             n_x = ((2+3*obj.lambda) / obj.lambda);
%             n_y = ((2+obj.lambda) / obj.lambda);

            % Derivative Wetting phase relative permeability 
            dkr(:,1) = (obj.n_w).* ((s-Phases(2).sr)/(1-Phases(2).sr)).^(obj.n_w-1).*(1/(1-Phases(2).sr));    
            dkr(s < Phases(2).sr, 1) = 0;
            
            % Derivative Non-Wetting phase relative permeability 
            dkr(:,2) = ((-2 + 2 * s)/((1-Phases(2).sr).^2)) .* (1 - ((s-Phases(2).sr) / (1 - Phases(2).sr)).^(obj.n_nw)) + ((1 - s)/(1 - Phases(2).sr)).^(2) .* -((obj.n_nw/(1-Phases(2).sr)) .* ((s-Phases(2).sr) / (1 - Phases(2).sr)).^(obj.n_nw-1));
            %dkr(s > 1 - Phases(1).sr, 2) = 0;
            
            %dkr(s < Phases(2).sr, 2) = 0;
            %dkr(s > 1 - Phases(1).sr, 1) = 0;
            
%             % lambda formulations
%             dkr(:,3) = (n_x).* ((s-Phases(1).sr)/(1-Phases(1).sr)).^(n_x-1).*(1/(1-Phases(1).sr));
%             dkr(:,4) = ((-2 + 2 * s)/((1-Phases(1).sr).^2)) .* (1 - ((s-Phases(1).sr) / (1 - Phases(1).sr)).^(n_y)) + ((1 - s)/(1 - Phases(1).sr)).^(2) .* -((n_y/(1-Phases(1).sr)) .* ((s-Phases(1).sr) / (1 - Phases(1).sr)).^(n_y-1));
%             dkr(s < Phases(1).sr, 3) = 0;
%             dkr(s > 1 - Phases(2).sr, 4) = 0;
%             dkr(s < Phases(1).sr, 4) = 0;
%             dkr(s > 1 - Phases(2).sr, 3) = 0;
             
        end
    end
end
   