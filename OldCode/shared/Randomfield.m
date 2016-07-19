% 2D random permeability field generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Tom Postma
%TU Delft
%Created: January 2016
%Last modified: January 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K] = Randomfield(Grid)
% Gaussian distribution function and autocovariance functions
    
    Nx      = Grid.Nx;
    Ny      = Grid.Ny;
    K_av    = Grid.Kav;
    L_cx    = Grid.Lcx;
    L_cy    = Grid.Lcy;
    V_dp    = Grid.Vdp;
       
    x   = linspace(-1/2,1/2,Nx);
    y   = linspace(-1/2,1/2,Ny);

    [X,Y] = meshgrid(x,y); 

    s   = -log(1-V_dp);
    mu  = log(K_av)-(s^2)/2;
    Z   = s*randn(Nx,Ny);
    
    % Gaussian filter
    a = X.^2/(L_cx^2/2);
    b = Y.^2/(L_cy^2/2) ; 
    F = exp(-(a+b));
    
    % correlated surface generation
    f   = 2.0/sqrt(pi)*1/sqrt(Nx*Ny)/sqrt(L_cx)/sqrt(L_cy)*ifft2(fft2(Z).*fft2(F));
    K   = exp(mu+real(f));
    
end




