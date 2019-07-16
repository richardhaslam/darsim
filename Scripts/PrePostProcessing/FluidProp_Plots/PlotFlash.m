% Script
%z = 0:0.01:1;
z = 0:0.01:1;
p = 0;
x1v = zeros(length(z), 1);
x1l = zeros(length(z), 1);
x2v = zeros(length(z), 1);
x2l = zeros(length(z), 1);
ni = zeros(length(z), 1);
for i=1:length(z)
    [x1v(i), x1l(i), x2v(i), x2l(i), ni(i)] = BOFlash(z(i), p);
end
figure(1)
plot(z, x1v);
xlabel('z1');
ylabel('x1v');
figure(2)
plot(z, x1l);
xlabel('z1');
ylabel('x1l');
figure(3)
plot(z, x2v);
xlabel('z1');
ylabel('x2v');
figure(4)
plot(z, x2l);
xlabel('z1');
ylabel('x2l');
figure(5)
plot(z, ni);
xlabel('z1');
ylabel('ni');

z = 0.5;
p = 0:0.01:1;
x1v = zeros(length(z), 1);
x1l = zeros(length(z), 1);
x2v = zeros(length(z), 1);
x2l = zeros(length(z), 1);
ni = zeros(length(z), 1);
for i=1:length(p)
    [x1v(i), x1l(i), x2v(i), x2l(i), ni(i)] = BOFlash(z, p(i));
end
figure(10)
plot(p, x1v);
xlabel('z1');
ylabel('x1v');
figure(20)
plot(p, x1l);
xlabel('z1');
ylabel('x1l');
figure(30)
plot(p, x2v);
xlabel('z1');
ylabel('x2v');
figure(40)
plot(p, x2l);
xlabel('z1');
ylabel('x2l');
figure(50)
plot(p, ni);
xlabel('z1');
ylabel('ni');
