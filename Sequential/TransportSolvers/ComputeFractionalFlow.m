%Fractional Flow Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fw, dFw] = ComputeFractionalFlow(s, Fluid, Grid, Utot)
%Compute Mobilities
[Mw, Mo] = Mobilities(s, Fluid);
Mt = Mw + Mo;
%Flux function and derivatives
fw = Mw./Mt;
dfw = Derivative(s,Fluid);
if (isempty(Fluid.Pc))
    % Without capillarity
    Fw = fw;
    dFw = dfw;
else
    % With capillarity
    Utot.x;
    Utot.y;
    fw = reshape(fw, Grid.Nx, Grid.Ny);    
    CapFx(2:Nx+1,:) = (C(:,:) - C(:,:))./Grid.dx;
    CapFy(:,2:Ny+1) = (C(:,:) - C(:,:))./Grid.dy;
    Fw(:,:) = fw(:,:) + CapFx + CapFy;
    Fw = reshape(Fw, Grid.N, 1);
end
end