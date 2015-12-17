function CheckEllipticSystem(A, Nx)
A = sign(A);
main = diag(A);
x1 = diag(A, 1);
x2 = diag(A, -1);
y1 = diag(A, Nx);
y2 = diag(A, -Nx);
[a, b] = find(main < 0);
[c, d] = find(x1 > 0);
[e, f] = find(x2 > 0);
[g, h] = find(y1 > 0);
[m, n]= find(y2 > 0);
end