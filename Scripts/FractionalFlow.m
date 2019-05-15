% f = (rho1 * Mob1) / (rho1 * Mob1 + rho2 * Mob2);
% df = (rho1 * dMob1 * (rho1*Mob1 + rho2*Mob2) - (rho1*dMob1 + rho2*dMob2) * rho1*Mob1) / (rho1 * Mob1 + rho2 * Mob2)^2
% Num = A - B;
% A = (rho1 * dMob1 * (rho1*Mob1 + rho2*Mob2)
% dA = rho1 * ddMob1 * (rho1*Mob1 + rho2*Mob2) + rho1 * dMob1 * (rho1*dMob1 + rho2 * dMob2)
% B = (rho1*dMob1 + rho2*dMob2) * rho1*Mob1
% dB = (rho1*ddMob1 + rho2 * ddMob2) *  rho1*Mob1 + (rho1*dMob1 + rho2*dMob2) * rho1*dMob1
% Den = (rho1 * Mob1 + rho2 * Mob2)^2;
% dDen = 2 * (rho1 * Mob1 + rho2 * Mob2) * (rho1 * dMob1 + rho2 * dMob2) 

s = [0:0.001:1]';
mu = [1e-3, 1.5e-3];
rho = [1000*ones(length(s), 1), 100*ones(length(s), 1)];
Phases(1).sr = 0.0;
Phases(2).sr = 0.0;
RelPerm = relperm_model_quadratic();
kr = RelPerm.ComputeRelPerm(Phases, s);
dkr = RelPerm.ComputeDerivative(Phases, s);
ddkr = RelPerm.ComputeSecondDerivative(Phases, s);

Mob = 0 * kr;
dMob = 0 * kr;
ddMob = 0 * kr;
for i = 1:2
    Mob(:,i) = kr(:, i) ./ mu(i);
    dMob(:,i) = dkr(:,i) ./ mu(i);
    ddMob(:,i) = ddkr(:,i) ./ mu(i);
end

f = rho(:, 1) .* Mob(:,1) ./ (rho(:, 1) .* Mob(:,1) + rho(:, 2) .* Mob(:, 2));
num = rho(:, 1) .* Mob(:, 1);
dnum = rho(:, 1) .* dMob(:,1);
den = sum(rho .* Mob, 2);
dden = sum(rho .* dMob, 2);
df = (dnum .* den - dden .* num) ./ den.^2;

A = rho(:, 1) .* dMob(:, 1) .* sum(rho .* Mob, 2);
dA = rho(:, 1) .* ddMob(:,1) .* sum(rho.*Mob, 2) + rho(:,1) .* dMob(:,1) .*  sum(rho.*dMob, 2);
B = sum(rho.*dMob, 2) .* rho(:, 1) .* Mob(:,1);
dB = sum(rho.*ddMob, 2) .*  rho(:, 1).*Mob(:, 1) + sum(rho.*dMob, 2) .* rho(:, 1) .* dMob(:, 1);
num = A - B;
dnum = dA - dB;
den = sum(rho.*Mob, 2).^2;
dden = 2 * sum(rho.*Mob, 2) .* sum(rho.*dMob, 2);
ddf = (dnum .* den - dden .* num) ./ den.^2;

figure(1)
plot(s, f);
title('f(s)');
figure(2)
plot(s, df);
title('df(s)');
figure(3)
plot(s, ddf);
title('ddf(s)');