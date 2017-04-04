K = load('../Permeability/spe10/Permx.dat');
K = K(4:end);

kspe10 = zeros(length(K),1);
for i=1:84
    kspe10((i-1)*220*60 + 1:220*60*(i+1)) = K((85-i-1)*220*60+1:220*60*(85-i+1));
end
kspe10(84*220*60 + 1:220*60*85) = K(1:220*60);
kspe10 = [60; 220; 85; kspe10];

save('../Permeability/spe10/Permxflip.dat', 'kspe10', '-ascii');