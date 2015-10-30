function Df=Derivative(s, Fluid)
%dfdS
%This function computes the derivative of the fractional flow function with
%respect to S at a given saturation value.
S = (s-Fluid.swc)/(1-Fluid.swc-Fluid.sor); % Rescale saturations
    swc=Fluid.swc;
    sor=Fluid.sor;
    muw=Fluid.muw;
    muo=Fluid.muo;
if (strcmp(Fluid.RelPerm,'Quadratic')==1)
    Df=-(1-swc-sor).^(-1).*((2*muw*muo.*(S-1).*S)./(muw*(S-1).^2+muo.*S.^2).^2);
else
    Df=(1-swc-sor).^(-1)*(muw*muo)./(muw.*(S-1)-muo.*S).^2;
end
end