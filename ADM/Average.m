function [P_, S_] = Average(P, S, CoarseGrid, FineGrid)
%Average Solution on Coarse Blocks
Nc = CoarseGrid.Nx*CoarseGrid.Ny;
N = FineGrid.Nx*FineGrid.Ny;
p = reshape(P, N, 1);
s = reshape(S, N, 1);
p_ = p;
s_ = s;
for c=1:Nc
    if CoarseGrid.Active(c) == 1
        Imin = CoarseGrid(1).I(c) - floor((CoarseGrid(1).CoarseFactor(1) - 1)/2);
        Imax = CoarseGrid(1).I(c) + ceil((CoarseGrid(1).CoarseFactor(1) - 1)/2);
        Jmin = CoarseGrid(1).J(c) - floor((CoarseGrid(1).CoarseFactor(2) - 1)/2);
        Jmax = CoarseGrid(1).J(c) + ceil((CoarseGrid(1).CoarseFactor(2) - 1)/2);
        i=Imin:Imax;
        j=Jmin:Jmax;
        [I, J] = meshgrid(i, j);
        pairs = [I(:), J(:)];
        z = pairs(:,1) + (pairs(:,2)-1)*FineGrid.Nx;
        p_av = mean(p(z));
        s_av = mean(s(z));
        p_(z) = p_av;
        s_(z) = s_av;
    end
end
P_ = reshape(p_, FineGrid.Nx, FineGrid.Ny);
S_ = reshape(s_, FineGrid.Nx, FineGrid.Ny);
end