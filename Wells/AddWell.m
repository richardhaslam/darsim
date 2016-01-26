function [A, q] = AddWell(A, q, Well, K)
a = Well.cells.cells;
for i=1:length(a)
    A(a(i),a(i)) = A(a(i),a(i))+ Well.PI*K(a(i));
    q(a(i)) = Well.PI*K(a(i))*Well.p;
end
end