function [x1v, x1l, x2v, x2l, ni] = BOFlash(z, p)
% 
rho1ST = 100;
rho2ST = 800;
Rs = 0.2*p + 0.2;
k1 = (rho1ST * Rs + rho2ST) ./ (Rs * rho1ST);
k2 = 0;
% Initial guess
x1v = 1;
x1l = 0.5;
x2v = 0;
x2l = 0.5;
ni = 0.5;

% Initial residual
R = zeros(5,1);
R(1) = x1v - k1 * x1l;
R(2) = x2v - k2 * x2l;
R(3) = z - (1-ni) * x1l - ni * x1v;
R(4) = 1 - z - (1-ni) * x2l - ni * x2v;
R(5) = x1v - x1l + x2v - x2l;

while norm(R) > 1e-6    
    A = zeros(5);
    % Equation 1
    % dF1/dx1v
    A(1, 1) = 1;
    % dF1/dx1l
    A(1, 2) = -k1;
    % Equation 2
    % dF2/dx2v
    A(2, 3) = 1;
    % dF2/dx2l
    A(2, 4) = -k2;
    % Equation 3
    % dF3/dx1v
    A(3, 1) = - ni;
    % dF3/dx1l
    A(3, 2) = -(1 - ni);
    % dF3dni
    A(3, 5) = x1l - x1v;
    % Equation 4
    % dF4/dx2v
    A(4, 3) = - ni;
    % dF4/dx2l
    A(4, 4) = - (1 - ni);
    % dF4dni
    A(4, 5) = x2l - x2v;
    % Equation 5
    % dF5/dx1v
    A(5, 1) = 1;
    % dF5/dx1l
    A(5, 2) = -1;
    % dF5/dx2v
    A(5, 3) = 1;
    % dF5/dx2l
    A(5, 4) = -1;
    
    Sol = - A\R;
    
    x1v = x1v + Sol(1);
    x1l = x1l + Sol(2);
    x2v = x2v + Sol(3);
    x2l = x2l + Sol(4);
    ni =  ni + Sol(5);
    % Update residual
    R(1) = x1v - k1 * x1l;
    R(2) = x2v - k2 * x2l;
    R(3) = z - (1-ni) * x1l - ni * x1v;
    R(4) = 1 - z - (1-ni) * x2l - ni * x2v;
    R(5) = x1v - x1l + x2v - x2l;    
end
    if (x1l > z)
        x1l = z;
        x2l = 1-z;
        x1v = 0;
        x2v = 0;
        ni = 0;
    end
end