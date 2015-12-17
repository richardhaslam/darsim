%Phase mobilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mw, Mo, dMw, dMo]=Mobilities(s, Fluid)
%Phase Mobilities
S = (s-Fluid.swc)/(1-Fluid.swc-Fluid.sor); % Rescale saturations
if (strcmp(Fluid.RelPerm,'Quadratic')==1)

    Mw = S.^2 /Fluid.muw;        % Water mobility
    Mo = (1-S).^2/Fluid.muo;     % Oil mobility
    dMw = (1-Fluid.swc-Fluid.sor)^(-1)*2*S/Fluid.muw;
    dMo = -(1-Fluid.swc-Fluid.sor)^(-1)*2*(1-S)/Fluid.muo;
	
elseif (strcmp(Fluid.RelPerm, 'Foam')==1)
    % w = light phase (gas)
    % o = heavy phase


    swc     = Fluid.swc;
    sor     = Fluid.sor;   
    krwe    = Fluid.krwe;
    kroe    = Fluid.kroe;
    nw      = Fluid.nw;
    no      = Fluid.no;
    muw     = Fluid.muw;
    muo     = Fluid.muo;
    epdry   = Fluid.epdry;
    fmmob   = Fluid.fmmob;
    fmdry   = Fluid.fmdry;


    Mx     = (krwe * (S.^nw))/muw;
    Mo     = (kroe * ((1-S).^no))/muo;

    FW     = 0.5 + (1/pi) .* atan(epdry.*((1-s)-fmdry));
    FI     = ((1+(fmmob.*FW)).^(-1));
    
    Mw     = Mx .* FI;
    
    d1     =    (1-swc-sor)^(-1) * krwe * nw .*   ((S).^(nw-1)) ./ muw;
    d2     =    (fmmob.*epdry) ./ (pi.*(epdry^2 .* (fmdry+s-1).^2 + 1) .* ((-1/pi) .* fmmob .* atan(epdry*(fmdry+s-1))+0.5.*fmmob+1).^2);

    dMw = d1 .* FI + Mx .* d2;
    dMo = -1 * (1-swc-sor)^(-1) * kroe * no .* ((1-S).^(no-1)) ./ muo;	
else
    Mw=S/Fluid.muw;              %Water mobility
    Mo=(1-S)/Fluid.muo;          %Oil mobility
    dMw = ones(length(s),1)*(1-Fluid.swc-Fluid.sor)^(-1)/Fluid.muw;
    dMo = -ones(length(s),1)*(1-Fluid.swc-Fluid.sor)^(-1)/Fluid.muo;
end
end


