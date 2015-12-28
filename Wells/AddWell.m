function [A, q] = AddWell(A, q, Well, K, Nx)
a=Well.x+(Well.y-1)*Nx;
if (Well.Dirichlet == 1)
    for i=1:length(a)
        A(a(i),:) = 0;
        A(a(i), a(i)) = 1; 
        q(a(i)) = Well.p;
    end
else
    for i=1:length(a)
        A(a(i),a(i))=A(a(i),a(i))+ Well.PI*K(1,Well.x(i),Well.y(i));
        q(a)=Well.PI*K(1,Well.x(i),Well.y(i))*Well.p;
    end
end
end