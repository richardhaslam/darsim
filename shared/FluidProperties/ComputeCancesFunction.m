%Compute Cances function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C] = ComputeCancesFunction(Fluid, S, K, Grid)
%s = 0.1:0.01:1;
[Mw, Mo] = Mobilities(S, Fluid);
Mt = Mw + Mo;
[~, dPc] = ComputePc(S, Fluid);
integrand = K.*(Mw.*Mo)./Mt.*dPc;
integrand = reshape(integrand, Grid.N, 1);
C = zeros(Grid.N, 1);
for i=1:size(Grid.N)
    C(i) = Integral(integrand(i), 0, S);
end
C = reshape(C, Grid.Nx, Grid.Ny);
%figure(177)
%plot(s, C);
%xlabel('water saturation');
%ylabel('Cances function');
end