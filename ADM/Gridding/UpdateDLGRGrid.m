%ADM - Update ADM grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NewDLGRGrid = UpdateDLGRGrid(DLGRGrid, FineGrid, Nc, Nf, Nx)
level = DLGRGrid.level(end);
[r, c] = size(DLGRGrid.Father);
%Add cells to DLGR grid
Field1 = 'N';  Value1 = DLGRGrid.N;
Field2 = 'I';  Value2 = zeros(Nf, 1);
Field3 = 'J';  Value3 = zeros(Nf, 1);
Field4 = 'CoarseFactor'; Value4 = zeros(Nf, 1); 
Field5 = 'CellIndex'; Value5 = zeros(Nf, 1);
Field6 = 'Father'; Value6 = zeros(Nf, c);
Field7='Centers'; Value7 = zeros(Nf,c);
NewDLGRGrid = struct(Field1,Value1,Field2,Value2,Field3,Value3,Field4,Value4,Field5,Value5, Field6, Value6, Field7, Value7);
NewDLGRGrid.I(1:Nx) = DLGRGrid.I(1:Nx);
NewDLGRGrid.J(1:Nx) = DLGRGrid.J(1:Nx);
NewDLGRGrid.CoarseFactor(1:Nx) = DLGRGrid.CoarseFactor(1:Nx);
NewDLGRGrid.CellIndex(1:Nx) = DLGRGrid.CellIndex(1:Nx);
NewDLGRGrid.Father(1:Nx, :) = DLGRGrid.Father(1:Nx,:);
NewDLGRGrid.Centers(1:Nx, :) = DLGRGrid.Centers(1:Nx,:);
NewDLGRGrid.level(1:Nx) = DLGRGrid.level(1:Nx);
NewDLGRGrid.N(level+1) = 0;
h=Nx+1;
for i=1:FineGrid.Nx*FineGrid.Ny
    if FineGrid.Active(i) == 0
        for j=Nx+1:Nc
            if FineGrid.Father(i, level) == DLGRGrid.CellIndex(j)
                NewDLGRGrid.I(h) = FineGrid.I(i);
                NewDLGRGrid.J(h) = FineGrid.J(i);
                NewDLGRGrid.CoarseFactor(h,1) = FineGrid.CoarseFactor(1);
                NewDLGRGrid.CoarseFactor(h,2) = FineGrid.CoarseFactor(2);
                NewDLGRGrid.CellIndex(h) = i;
                NewDLGRGrid.level(h) = level - 1;
                NewDLGRGrid.Father(h,:) = FineGrid.Father(i,:);
                NewDLGRGrid.Centers(h,:) = FineGrid.Centers(i,:);
                NewDLGRGrid.N(level) = NewDLGRGrid.N(level) + 1;
                h = h + 1;
            end
        end
    end
end
end