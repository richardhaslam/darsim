%Plotting - fractional flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FractionalFlowPlot(Fluid)
%Plot fractional flow function and its derivatives
X = linspace(Fluid.swc,1 - Fluid.sor,100);
[Mw, Mo, dMw, dMo]=Mobilities(X, Fluid);
fw=Mw./ (Mw+Mo);
df=Derivative(X, Fluid);
Ddf=Derivative_2nd(X, Fluid);

%%%Fractional flow
figure(1001);
plot(1-X, fw);
xlabel('S');
ylabel('fw');
figure(1002);
plot(1-X, df);
xlabel('S');
ylabel('dfw');
figure(1003);
plot(1-X, Ddf);
xlabel('S');
ylabel('Ddfw');

%%%Mobilities
figure(1004)
plot(X, Mw);
xlabel('S');
ylabel('Mw');
figure(1005)
plot(X, dMw);
xlabel('S');
ylabel('dMw');
figure(1006)
plot(1-X, Mo);
xlabel('S');
ylabel('Mo');
figure(1007)
plot(1-X, dMo);
xlabel('S');
ylabel('dMo');
end

