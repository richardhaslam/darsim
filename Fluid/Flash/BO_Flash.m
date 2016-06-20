function [x] = BO_Flash(p)

%% - Description
PwfValue = 10e6;                                 %HARD CODED to set up Pdim

%% - Dimensionless pressure!
Pdim = p/max(PwfValue);                  %Find dimensionless pressure to perform flash      

%% - Solve for x's
x(:,2) = 1 - (800./(800 + 100*(0.2*Pdim(:,1) + 0.2)));  %More pressure less comp 1 (oil) in phase 1 (oil) (actually more gas pushed in really)
x(:,1) = 1 - 0;                              %Gas is all gas

end

