%%% Test Rachford-Rice convergence
z = zeros(1, 2);
z(1) = 0.7499;
z(2) = 1 - z(1);
k = [1.5, 0.5];

[fv, x, y] = RachfordRice(z, k);
