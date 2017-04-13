%% J function Pc curves
s = 0:0.01:1;
swr = 0.1;
snwr = 0.1;
alpha  = 4.361e-2;
K = 1e-13;
por = 0.2;

% Rescale sat
S = (s - snwr)/(1 - snwr) + 0.1;
S(s<snwr) = 0.1;

% J function
J = 0.1.*(S).^(-0.5);
dJ = - 0.1*0.5*(S).^(-1.5);

%Pc
Pc = alpha * (por./K).^(0.5) * J;
dPc = alpha .* (por./K).^(0.5) .* dJ;

figure(1)
title('Pc vs Sw')
plot(s, Pc, 'red');
xlim([0 1]);

figure(2)
title('dPc vs Sw')
plot(s, dPc, 'blue');
xlim([0 1]);

