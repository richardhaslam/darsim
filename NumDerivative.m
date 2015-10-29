function Df=NumDerivative(S, Fluid)
%dfdS
%This function computes the derivative of the fractional flow function with
%respect to S at a given saturation value.
Delta=0.0001;
if (S>(1-Delta))
    S2=S;
    S1=S-Delta;
else
    S1=S;
    S2=S+Delta;
end
[Mw,Mo]=Mobilities(transpose([S1, S2]),Fluid); % compute mobilities
fw = Mw./(Mw+Mo);            % compute fractional flow function
Df=(fw(2)-fw(1))/Delta;
end