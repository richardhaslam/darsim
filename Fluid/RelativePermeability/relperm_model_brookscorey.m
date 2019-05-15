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
      n_w = 2;
      n_nw = 2; 
        %lambda = 3;

    end
        
    methods
        
        function kr = ComputeRelPerm(obj, Phases, s)     
            % phase 1 is non-wetting here               
            S = (s-Phases(2).sr)/(1-Phases(2).sr);
            S = max(S, 0);
            
            % wetting phase rel perm
            kr(:,2) = (S).^(obj.n_w);           
            %S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            %kr(s < Phases(2).sr, 1) = 0;
           
            %Non-Wetting phase relative permeability
            kr(:,1) = (1 - S).^ 2 .*( 1 - (S).^(obj.n_nw));
            
            %kr(:,2) = ( (1-s) / (1-Phases(2).sr-Phases(1).sr) ) .^ 2 .* ( 1 - ((s-Phases(2).sr)/(1-Phases(2).sr-Phases(1).sr)).^(obj.n_nw));
            %kr(:,2) = (1- S) .^ 2  -  ((1- S).^ 2).*((S).^(obj.n_nw));
            %kr(s < Phases(1).sr, 1) = 1;
            %kr(s < Phases(2).sr, 1) = 1
            %kr(s > 1 - Phases(1).sr, 1) = 0;
            
            kr(s < Phases(2).sr, 1) = 1;    %DISCUSS
            kr(s > 1- Phases(1).sr, 1) = 0;
            kr(s > 1- Phases(1).sr, 2) = 1;               
            kr(s < Phases(2).sr, 2) = 0;
            
            %kr(:,1) = smooth(kr(:,1));
            
            %kr(s > 1 - Phases(1).sr, 2) = 0;
            %kr(s > 1 - Phases(1).sr,1) = 1;
            
        end
        function dkr = ComputeDerivative(obj, Phases, s)
            %S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            %S = (s-Phases(1).sr)/(1-Phases(1).sr);
            S = (s-Phases(2).sr)/(1-Phases(2).sr);
            
            % Derivative Wetting phase relative permeability 
            dkr(:,2) = ( obj.n_w .* (S) .^ (obj.n_w-1) ) ./ (1-Phases(2).sr);

            %dkr(s < Phases(1).sr, 1) = 0;
            %dkr(:,1) = ( obj.n_w .* (S) .^ (obj.n_w-1) ) ./ (1-Phases(2).sr-Phases(1).sr);    
            %dkr(:,1) = ( obj.n_w .* (S) .^ (obj.n_w-1) ) ./ (1-Phases(1).sr);
            
            % Derivative Non-Wetting phase relative permeability 
            dkr(:,1) = ((-2 + 2 * s)/((1-Phases(2).sr).^2)) .* (1 - (S).^(obj.n_nw)) + ((1 - s)/(1 -Phases(2).sr)).^(2) .* -((obj.n_nw/(1-Phases(2).sr)) .* (S).^(obj.n_nw-1));
            
            dkr(s < Phases(2).sr, 1) = 0; %DISCUSS
            dkr(s > 1- Phases(1).sr, 1) = 0;
            dkr(s > 1- Phases(1).sr, 2) = 0;
            dkr(s < Phases(2).sr, 2) = 0;
            
            %dkr(:,2) = ((-2 + 2 * s)/((1-Phases(1).sr-Phases(1).sr).^2)) .* (1 - (S).^(obj.n_nw)) + ((1 - s)/(1 - Phases(2).sr-Phases(1).sr)).^(2) .* -((obj.n_nw/(1-Phases(2).sr-Phases(1).sr)) .* (S).^(obj.n_nw-1));
            %dkr(:,2) = ((-2 + 2 * s)/((1-Phases(1).sr).^2)) .* (1 - (S).^(obj.n_nw)) + ((1 - s)/(1 -Phases(1).sr)).^(2) .* -((obj.n_nw/(1-Phases(1).sr)) .* (S).^(obj.n_nw-1));
            
            %dkr(s > 1 - Phases(1).sr, 2) = 0;
            %kr(s > 1 - Phases(1).sr, 2) = 0;
            %kr(s > 1 - Phases(1).sr,1) = 1;
            %dkr(imag(dkr) ~= 0) = 0;              % set imaginary numbers to zero in derivative
            %dkr(s < Phases(2).sr, 2) = 0;
            %dkr(s > 1 - Phases(1).sr, 1) = 0;
                       
        end
    end
end
   