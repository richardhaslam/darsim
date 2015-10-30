function PI = ComputeProductivityIndex(r_w, Kx, Ky, Dx, Dy, h)
%Computes productivity index of a well
%1. Compute equivalent well block radius
Num =((Ky/Kx)^0.5*Dx^2+(Kx/Ky)^0.5*Dy^2)^0.5;
Den = (Ky/Kx)^0.25 + (Kx/Ky)^0.25;
r_o = 0.28*Num/Den;
%2. Compute productivity index
PI = 2*pi*h/log(r_o/r_w);
end