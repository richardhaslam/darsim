% Let s check these dependencies 
index = 1;
z1 = [0:0.1:1]';
Statu.x = [z1 0 0 z1 0 z1];
Status.ni = 0.5*ones(length(z1),1);

for i=1:length(z1)
    s
    Physical = sum(Status.z, 2);
    if Physical(i) == 1
        [Status SinglePhase] = SimpleFlashCalculator(Status);
    else
        Status.ni(i) = -1;
        Status.x(i, :) = [-1 -1 -1 -1 -1 -1];
    end
end
