function [x] = BO_Flash(p, Fluid)

%% - Description

%% - Dimensionless pressure!
Pdim = p/Fluid.Pref;                  %Find dimensionless pressure to perform flash      

%% - Solve for x's
x(:,2) = 1 - (800./(800 + 100*(0.2*Pdim(:,1) + 0.2)));  %More pressure less comp 1 (oil) in phase 1 (oil) (actually more gas pushed in really)
x(:,1) = 1 - 0;                              %Gas is all gas

end

