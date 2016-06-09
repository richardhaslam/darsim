%Phase mobilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mnw, Mw, dMnw, dMw, app]=Mobilities(s, Fluid)

%Phase Mobilities
S = (s-Fluid.sr(2))/(1-Fluid.sr(2)-Fluid.sr(1)); % Rescale saturations

if (strcmp(Fluid.RelPerm,'Quadratic')==1)

    Mnw = S.^2 /Fluid.mu(2);        % Non-wetting mobility
    Mw = (1-S).^2/Fluid.mu(1);     % Wetting mobility
    dMnw = (1-Fluid.sr(2)-Fluid.sr(1))^(-1)*2*S/Fluid.mu(2);
    dMw = -(1-Fluid.sr(2)-Fluid.sr(1))^(-1)*2*(1-S)/Fluid.mu(1);
    
elseif (strcmp(Fluid.RelPerm, 'Linear')==1)
    Mnw=(S/Fluid.mu(2));                %Nonwetting mobility
    Mw=((1-S)/Fluid.mu(1));             %Wetting mobility
    dMnw = (ones(size(Mnw))).*(1-Fluid.sr(2)-Fluid.sr(1)).^(-1)./Fluid.mu(2);
    dMw = -(ones(size(Mw))).*(1-Fluid.sr(2)-Fluid.sr(1)).^(-1)./Fluid.mu(1);

elseif (strcmp(Fluid.RelPerm, 'Corey')==1)
    Mnw = Fluid.kre(2)*S.^(Fluid.n(2))/Fluid.mu(2);
    Mw = Fluid.kre(1)*(1-S).^(Fluid.n(1))/Fluid.mu(1);
    dMnw = Fluid.n(2)*(1/(1-Fluid.sr(2)-Fluid.sr(1)))*Fluid.kre(2)*S.^(Fluid.n(2)-1)/Fluid.mu(2);
    dMw = Fluid.n(1)*(1/(1-Fluid.sr(2)-Fluid.sr(1)))*Fluid.kre(1)*(1-S).^(Fluid.n(1)-1)/Fluid.mu(1);	
    
elseif (strcmp(Fluid.RelPerm, 'Foam')==1)
    % w = light phase (gas)
    % o = heavy phase
    epdry   = Fluid.epdry;
    fmmob   = Fluid.fmmob;
    fmdry   = Fluid.fmdry;


    Mx     = (Fluid.kre(2) * (S.^Fluid.n(2)))/Fluid.mu(2);
    Mw     = (Fluid.kre(1) * ((1-S).^Fluid.n(1)))/Fluid.mu(1);

    FW     = 0.5 + (1/pi) .* atan(epdry.*((1-s)-fmdry));
    FI     = ((1+(fmmob.*FW)).^(-1));
    
    Mnw     = Mx .* FI;
    
    app             = Fluid.mu(2) ./ FI;
    app((round(app,2)==0.5)) = 0;
    
    
    d1     =    (1-Fluid.sr(2)-Fluid.sr(1))^(-1) * Fluid.kre(2) * Fluid.n(2) .*   ((S).^(Fluid.n(2)-1)) ./ Fluid.mu(2);
    d2     =    (fmmob.*epdry) ./ (pi.*(epdry^2 .* (fmdry+s-1).^2 + 1) .* ((-1/pi) .* fmmob .* atan(epdry*(fmdry+s-1))+0.5.*fmmob+1).^2);

    dMnw = d1 .* FI + Mx .* d2;
    dMw = -1 * (1-Fluid.sr(2)-Fluid.sr(1))^(-1) * Fluid.kre(1) * Fluid.n(1) .* ((1-S).^(Fluid.n(1)-1)) ./ Fluid.mu(1);	
end
end


