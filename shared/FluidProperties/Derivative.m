%Fractional flow derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Df = Derivative(s, Fluid)
%dfdS
%   This function computes the derivative of the fractional flow function with
%   respect to S at a given saturation value.
S = (s-Fluid.swc)/(1-Fluid.swc-Fluid.sor); % Rescale saturations
swc = Fluid.swc;
sor = Fluid.sor;
muw = Fluid.muw;
muo = Fluid.muo;

if (strcmp(Fluid.RelPerm,'Quadratic') == 1)

    Df = -(1-swc-sor).^(-1).*((2*muw*muo.*(S-1).*S)./(muw*(S-1).^2+muo.*S.^2).^2);
	
elseif (strcmp(Fluid.RelPerm, 'Foam') == 1)

    krwe    = Fluid.krwe;
    kroe    = Fluid.kroe;
    nw      = Fluid.nw;
    no      = Fluid.no;
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

    Df        = (dMw .* (Mw+Mo)...
                - (dMw+dMo).* Mw)...
                ./(Mw+Mo).^2;     
    
	
else
    Df = (1-swc-sor).^(-1)*(muw*muo)./(muw.*(S-1)-muo.*S).^2;
end
end