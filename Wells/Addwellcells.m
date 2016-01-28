%Perforated cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Well] = Addwellcells(Nx, Direction, I, Y, final)
if strcmp(Direction, 'vertical')
    Well.cells = I + (Y-1)*Nx;
else
    Yindexes = Y:1:final;
    Well.cells = ones(1, length(Yindexes))*I + (Yindexes - 1)*Nx;
end
end