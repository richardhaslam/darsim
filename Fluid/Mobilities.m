%Phase mobilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mw, Mnw, dMw, dMnw, app] = Mobilities(s, Fluid)

%Phase Mobilities
S = (s-Fluid.sr(1))/(1-Fluid.sr(1)-Fluid.sr(2)); % Rescale saturations

switch(Fluid.RelPerm)
    case('Quadratic')
        
        Mw = S.^2 /Fluid.mu(1);        % Wetting phase mobility
        Mnw = (1-S).^2/Fluid.mu(2);     % Non-Wetting phase mobility
        dMw = (1-Fluid.sr(2)-Fluid.sr(1))^(-1)*2*S/Fluid.mu(1);
        dMnw = -(1-Fluid.sr(2)-Fluid.sr(1))^(-1)*2*(1-S)/Fluid.mu(2);
        
    case('Linear')
        Mw = (S/Fluid.mu(1));                %Nonwetting mobility
        Mnw = ((1-S)/Fluid.mu(2));             %Wetting mobility
        dMw = (ones(size(Mnw))).*(1-Fluid.sr(2)-Fluid.sr(1)).^(-1)./Fluid.mu(1);
        dMnw = -(ones(size(Mw))).*(1-Fluid.sr(2)-Fluid.sr(1)).^(-1)./Fluid.mu(2);
        
    case('Corey')
        Mw = Fluid.kre(1)*S.^(Fluid.n(1))/Fluid.mu(1);
        Mnw = Fluid.kre(2)*(1-S).^(Fluid.n(2))/Fluid.mu(2);
        dMw = Fluid.n(1)*(1/(1-Fluid.sr(2)-Fluid.sr(1)))*Fluid.kre(1)*S.^(Fluid.n(1)-1)/Fluid.mu(1);
        dMnw = Fluid.n(2)*(1/(1-Fluid.sr(2)-Fluid.sr(1)))*Fluid.kre(2)*(1-S).^(Fluid.n(2)-1)/Fluid.mu(2);
        
    case('Foam')
        % w = light phase (gas)
        % o = heavy phase
        epdry   = Fluid.epdry;
        fmmob   = Fluid.fmmob;
        fmdry   = Fluid.fmdry;
        
        
        Mx     = (Fluid.kre(2) * (S.^Fluid.n(2)))/Fluid.mu(2);
        Mnw     = (Fluid.kre(1) * ((1-S).^Fluid.n(1)))/Fluid.mu(1);
        
        FW     = 0.5 + (1/pi) .* atan(epdry.*((1-s)-fmdry));
        FI     = ((1+(fmmob.*FW)).^(-1));
        
        Mw     = Mx .* FI;
        
        app             = Fluid.mu(2) ./ FI;
        app((round(app,2)==0.5)) = 0;
        
        
        d1     =    (1-Fluid.sr(2)-Fluid.sr(1))^(-1) * Fluid.kre(2) * Fluid.n(2) .*   ((S).^(Fluid.n(2)-1)) ./ Fluid.mu(2);
        d2     =    (fmmob.*epdry) ./ (pi.*(epdry^2 .* (fmdry+s-1).^2 + 1) .* ((-1/pi) .* fmmob .* atan(epdry*(fmdry+s-1))+0.5.*fmmob+1).^2);
        
        dMw = d1 .* FI + Mx .* d2;
        dMnw = -1 * (1-Fluid.sr(2)-Fluid.sr(1))^(-1) * Fluid.kre(1) * Fluid.n(1) .* ((1-S).^(Fluid.n(1)-1)) ./ Fluid.mu(1);
end
end


