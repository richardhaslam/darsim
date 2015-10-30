%Plotting - fractional flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FractionalFlowPlot(Fluid)
%Plot fractional flow function and its derivatives
Fluid.RelPerm='Quadratic';
X=linspace(0,1,100);
[Mw, Mo]=Mobilities(X, Fluid);
fw=Mw./ (Mw+Mo);
df=Derivative(X, Fluid);
Ddf=Derivative_2nd(X, Fluid);
figure(1);
plot(X, fw);
xlabel('S');
ylabel('fw');
figure(2);
plot(X, df);
xlabel('S');
ylabel('dfw');
figure(3);
plot(X, Ddf);
xlabel('S');
ylabel('Ddfw');
end

