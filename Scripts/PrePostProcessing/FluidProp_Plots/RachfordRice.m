% Solve Rachford-Rice equation
function [fv, x, y]= RachfordRice(z, k)
fv = 0.992;
converged = 0;
alpha = 1;
iterations = 0;
while ~converged
    %Finds hi for each component
    hi = (z .* k) ./ (fv .* (k - 1) + 1);
    %Finds the derivative of hi for each component
    dhi = (z .* (k - 1).^2) ./ ((fv .* (k - 1) + 1).^2);
    h = sum(hi, 2) - 1;
    dh = - sum(dhi, 2);
    
    %Update fv
    fvnew = alpha .* (-h ./ dh) + fv;
    
    fv = fvnew;
    if norm(h, inf) < 1e-10
        converged = 1;
    end
    iterations = iterations + 1;
end
disp(['converged in ', num2str(iterations)])
x = zeros(1, 2);
y = zeros(1, 2);
% comp 1 in liquid
x(1) = z(1) ./ (fv .* (k(1) - 1) + 1);
x(2) = 1 - x(1);
% comp 1 in vapour
y(1) = k(1) .* x(1);
y(2) = 1 - y(1);
end