%%% Test Rachford-Rice convergence
z = zeros(1, 2);
z(1) = 0.1;
z(2) = 1 - z(1);
k = [14.2890, 0.0];

[fv_sol, x, y] = RachfordRice(z, k);

% Plot h and dh as functions of fv
fv = 0:0.00001:1;
hi = zeros(length(fv), 2);
dhi = zeros(length(fv), 2);
for i=1:2
    hi(:,i) = (z(i) .* k(i)) *1e6 ./ ((fv .* (k(i) - 1) + 1)*1e6);
    %Finds the derivative of hi for each component
    dhi(:,i) = (z(i) .* (k(i) - 1).^2) *1e6 ./ (((fv .* (k(i) - 1) + 1).^2) *1e6);
end
h = (sum(hi, 2) - 1)*1e6;
dh = (- sum(dhi, 2))*1e6;
figure(1)
plot(fv, h);
title('h(fv)');
xlabel('fv');
ylabel('h');
figure(2)
plot(fv, dh);
title('dh(fv)');
xlabel('fv');
ylabel('dh');
