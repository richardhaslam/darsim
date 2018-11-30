% Average Permeability

Perm = load('../Permeability/98deg/Perm_98deg_1.txt');
K = Perm(4:end);
K = reshape(K, Perm(1), Perm(2), Perm(3));

N = [99, 99];
K = reshape(K(1:N(1), 1:N(2)), N(1), N(2));
K = reshape(K, prod(N), 1);

cf = [9, 9];
Nc = N./cf;

R = sparse(prod(Nc), prod(N));

for i=1:Nc(1)
    for j=1:Nc(2)
        indexC = i + Nc(1)*(j-1);
        Imin = (i-1)*cf(1)+1;
        Imax = i*cf(1);
        Jmin = (j-1)*cf(2)+1;
        Jmax = j*cf(2);
        If = Imin:Imax;
        Jf = Jmin:Jmax;
        [p, q] = meshgrid(If, Jf);
        pairs = [p(:), q(:)];
        %indexes of the fine cells
        IndexesF = pairs(:,1) + (pairs(:,2)-1)*N(1);
        R(indexC, IndexesF) = 1;
    end
end
K_rest = R * K;
Kav = R' * (K_rest ./ sum(R, 2));

K_coarse = R*Kav/prod(cf);

figure(1)
pcolor(reshape(log(K), N(1), N(2)));
figure(2)
pcolor(reshape(log(K_coarse), Nc(1), Nc(2)));

K_coarse = [Nc';1;K_coarse];

save('K_c.txt','K_coarse', '-ascii');