%Returns Permeability Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K]=ReadPerm(Grid, Problem)
%Returns Permeability Field
Nx=Grid.Nx;
Ny=Grid.Ny;
N=Nx*Ny;
switch (Problem)
    case ('SPE10T')
        file  = '../Permeability/kk_t.txt';      % HERE YOU MENTION THE NAME OF THE FILE
        field = load(file);     % HERE YOU LOAD IT
        fieldX = field(2:2:size(field));      % HEAR YOU JUST READ THE KX ONES
        fieldX = reshape(fieldX,[60 220])';     % HERE YOU MAKE IT 60x220
        Kx = reshape(fieldX(1:Nx,1:Ny)*10^(-12), N, 1);  % HERE YOU TAKE THE SIZE YOU LIKE
        Ky =  Kx;
        K=reshape([Kx, Ky]', 2, Grid.Nx, Grid.Ny);
    case('SPE10B')
        file  = '../Permeability/kk_b.txt';      % HERE YOU MENTION THE NAME OF THE FILE
        field = load(file);     % HERE YOU LOAD IT
        fieldX = field(2:2:size(field));      % HEAR YOU JUST READ THE KX ONES
        fieldX = reshape(fieldX,[60 220])';     % HERE YOU MAKE IT 60x220
        Kx = reshape(fieldX(1:Nx,1:Ny)*10^(-12), N, 1);  % HERE YOU TAKE THE SIZE YOU LIKE
        Ky =  Kx;
        K=reshape([Kx, Ky]', 2, Grid.Nx, Grid.Ny);
    case ('Homogeneous')
        Kx = ones(N,1)*10^(-12);
        Ky = ones(N,1)*10^(-12);
        K=reshape([Kx, Ky]', 2, Grid.Nx, Grid.Ny);
    case('DARSim1')
        file  = '../Permeability/K_D.txt';      % HERE YOU MENTION THE NAME OF THE FILE
        field = load(file);     % HERE YOU LOAD IT
        fieldX = field(1:size(field));      % HEAR YOU JUST READ THE KX ONES
        fieldX = reshape(fieldX,[54 216])';     % HERE YOU MAKE IT 60x220
        Kx = reshape(fieldX(1:Nx,1:Ny), N, 1);  % HERE YOU TAKE THE SIZE YOU LIKE
        Ky =  Kx;
        K=reshape([Kx, Ky]', 2, Grid.Nx, Grid.Ny);
    case('DARSim2')
        file  = '../Permeability/K_darsim_2.txt';      % HERE YOU MENTION THE NAME OF THE FILE
        field = load(file);     % HERE YOU LOAD IT
        fieldX = field(1:size(field));      % HEAR YOU JUST READ THE KX ONES
        fieldX = reshape(fieldX,[54 216])';     % HERE YOU MAKE IT 60x220
        Kx = reshape(fieldX(1:Nx,1:Ny), N, 1);  % HERE YOU TAKE THE SIZE YOU LIKE
        Ky =  Kx;
        K=reshape([Kx, Ky]', 2, Grid.Nx, Grid.Ny);
end
end