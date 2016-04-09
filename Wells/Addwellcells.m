%Perforated cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Well] = Addwellcells(Nx, Direction, I, Y, final)
global Errors
switch Direction
    case ('vertical')
        Well.cells = I + (Y-1)*Nx;
    case ('horizontal')
        Yindexes = Y:1:final;
        Well.cells = ones(1, length(Yindexes))*I + (Yindexes - 1)*Nx;
    otherwise
        disp('ERROR: unknown well direction!! Please select either vertical or horizontal');
        Errors = 1;
        Well = 'error';
end
end