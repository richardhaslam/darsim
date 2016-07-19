%Computes second derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ddf=Derivative_2nd(s, Fluid)
S = (s-Fluid.swc)/(1-Fluid.swc-Fluid.sor); % Rescale saturations
swc = Fluid.swc;
sor = Fluid.sor;
muw = Fluid.muw;
muo = Fluid.muo;
K = (1-swc-sor).^(-1);

switch Fluid.RelPerm 
    case('Quadratic')
        Den = (muw*(S-1).^2+muo*S.^2).^3;
        Num = 2*muw*muo*(muw*(2*S+1).*(S-1).^2+muo*S.^2.*(2*S-3));
        Ddf = K.*Num./Den;
    case('Linear')
        Ddf = -K*2*muw*muo*(muw-muo)./(muw.*(S-1)-muo.*S).^3;
end
end