%% %%%%%%% GRAVITATIONAL FORCES %%%%%%%%
% This script plots fractional flow curves of wetting and non-wetting
% phases in presence of viscous and buoyancy forces. Quadratic relative
% permeability and constant velocities are considered.

%% Reservoir properties
K = 1e-15;  % m^2
H = 1000;   % m
deltaP = 10e6;  % Pa

%% Phase properties
% Phase 1 = wetting phase, Phase 2 = non-wetting phase
mu = [1e-3, 1e-3];              % Pa * s
rho =[995, 1000];               % kg/m^3
DeltaRho = rho(1) - rho(2);     % kg/m^3

%Total velocity (for now equal to U_c)
U_t = K*deltaP/H;   % m/s
g = 9.81;           % m/s^2

s = 0:0.01:1;
%% Compute viscous fractional flow term
kr = zeros(2, length(s));
kr (1,:) = s.^2;
kr (2,:) = (1-s).^2;
M(1,:) = kr(1,:) / mu(1);
M(2,:) = kr(2,:) / mu(2);
Mt = M(1,:) + M(2,:);
V(1,:) = M(1,:)./ Mt;
V(2,:) = M(2,:)./ Mt;

%% Compute gravity fractional flow part
m = mu(2)/mu(1);
% Gravity number
C_g = K * g * DeltaRho/ (mu(1)*U_t);
G = M(1,:).* kr(2,:)./Mt;

f(1,:) = V(1,:) - G * C_g;
f(2, :) = V(2,:) + G * C_g;

disp(['m = mu_n/mu_w: ' num2str(m)]);
disp(['mC_g: ' num2str(m*C_g)]);

%% Plotting

%Viscous term
figure(1)
suptitle('Viscous term');
subplot(1,2,1)
plot(s, V(1,:), 'blue');
title('wetting phase');
ylabel('V_w');
xlabel('S_w');
%axis([0 1 -1.5 1.5]);
subplot(1,2,2)
plot(1-s, V(2,:), 'blue');
title('non-wetting phase');
ylabel('V_n');
xlabel('S_n');
%axis([0 1 -1.5 1.5]);
hold on

%Gravity term
figure(2)
suptitle('Gravity term');
subplot(1,2,1)
plot(s, G, 'green');
title('wetting phase');
ylabel('G_w');
xlabel('S_w');
%axis([0 1 -1.5 1.5]);
subplot(1,2,2)
plot(1-s, G, 'green');
title('non-wetting phase');
ylabel('G_n');
xlabel('S_n');
%axis([0 1 -1.5 1.5]);
hold on

%Total term
figure(3)
%suptitle('Total term');
%subplot(1,2,1)
plot(s, V(1,:), 'blue', 'LineWidth', 3);
hold on
plot(s, G, 'green', 'LineWidth', 3);
hold on
plot(s, f(1,:), 'red', 'LineWidth', 3);
title('wetting phase');
ylabel('f_w');
xlabel('S_w');
legend('v','g', 'f');
set(gca,'fontsize',24);
%axis([0 1 -1.5 1.5]);
figure(4)
%subplot(1,2,2)
plot(1-s, V(2,:), 'blue');
hold on
plot(1-s, G, 'green');
hold on
plot(1-s, f(2,:), 'red');
title('non-wetting phase');
ylabel('f_n');
xlabel('S_n');
legend('v','g', 'f');
set(gca,'fontsize',24);
%axis([0 1 -1.5 1.5]);

