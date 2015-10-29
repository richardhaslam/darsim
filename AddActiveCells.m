function [DLGRGrid, Tot] = AddActiveCells(DLGRGrid, Grid, Num, level)
N = Grid.Nx*Grid.Ny;
count = 0;
for i=1:N
    if(Grid.Active(i) == 1)
        h = Num + count + 1;
        DLGRGrid.I(h) = Grid.I(i);
        DLGRGrid.J(h) = Grid.J(i);
        DLGRGrid.CoarseFactor(h, 1) = Grid.CoarseFactor(1);
        DLGRGrid.CoarseFactor(h, 2) = Grid.CoarseFactor(2);
        DLGRGrid.CellIndex(h) = i;
        DLGRGrid.level(h) = level;
        DLGRGrid.Father(h,:) = Grid.Father(i,:);
        DLGRGrid.Centers(h,:) = Grid.Centers(i,:);
        count = count +1;
    end
end
Tot = Num + count;
end