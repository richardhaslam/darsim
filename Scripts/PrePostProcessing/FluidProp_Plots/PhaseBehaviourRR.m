%%% Test Rachford-Rice convergence
z = zeros(1, 2);
z(1) = 0.7498  +   7.8128e-05;
z(2) = 1 - z(1);
k = [1.5, 0.5];

[fv, x, y] = RachfordRice(z, k);
