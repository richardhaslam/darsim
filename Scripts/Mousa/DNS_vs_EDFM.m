close all; clear all; clc;
Directory = 'Input/';
Nx = [225;045;025;015];
Ny = [225;045;025;015];
Nz = [001;001;001;001];
Lx = 9.0; dx = Lx./Nx;
Ly = 9.0; dy = Ly./Ny;
Cfx = Nx(1)./Nx(2:end);
Cfy = Ny(1)./Ny(2:end);
Cfz = Nz(1)./Nz(2:end); Cfz(Nz(2:end)==1)=0;
Nt = 235;
for n = 1:length(Nx)
    xcm{n} = linspace(dx(n)/2 , Lx-dx(n)/2 , Nx(n))';
    ycm{n} = linspace(dy(n)/2 , Ly-dy(n)/2 , Ny(n))';
end

Error_P = zeros(Nt,length(Nx)-1);
Error_S = zeros(Nt,length(Nx)-1);

dt = 1e-4; %[day]
%% Looping over time
for t = 1:Nt
    P = cell(length(Nx),1);
    S = cell(length(Nx),1);
    
    fprintf('Checking timestep: %f\n', t);
    % Reading Files
    for n = 1:length(Nx)
        if n==1
            File = strcat(Directory, 'ImmHomo_DNS_' ,num2str(Nx(n)),'_Sol',num2str(t),'.txt');
        else
            File = strcat(Directory, 'ImmHomo_EDFM_',num2str(Nx(n)),'_Sol',num2str(t),'.txt');
        end
        FileID = fopen(File, 'r');
        Matrix = textscan(FileID, '%s');
        Matrix = str2double(Matrix{1,1});
        Matrix = reshape(Matrix,[3,Nx(n)*Ny(n)*Nz(n)])';
        fclose(FileID);
        P{n} = Matrix(:,2);
        S{n} = Matrix(:,3);
    end
    
    % Comparing with Reference
    P_ref = cell(length(Nx),1);
    S_ref = cell(length(Nx),1);
    for n = 2:length(Nx)
        P_ref{n} = zeros(Nx(n)*Ny(n)*Nz(n),1);
        S_ref{n} = zeros(Nx(n)*Ny(n)*Nz(n),1);
        for k = 1:Nz(n)
            for j = 1:Ny(n)
                for i =1:Nx(n)
                    I_EDFM = Nx(n)*Ny(n)*(k-1) + Nx(n)*(j-1) + i;
                    I_ref  = Nx(1)*Ny(1)*(Cfz(1)*(k-1)+round(Cfz(1)/2)) + Nx(1)*(Cfy(1)*(j-1)+round(Cfy(1)/2)) + Cfx(1)*(i-1)+round(Cfx(1)/2);
                    P_ref{n}(I_EDFM) = P{1}(I_ref);
                    S_ref{n}(I_EDFM) = S{1}(I_ref);
                end
            end
        end
%         Error_P(t,n-1) = norm( (P_ref{n} - P{n}) ) / norm( P_ref{n} ) / length(P_ref{n});
%         Error_S(t,n-1) = norm( (S_ref{n} - S{n}) ) / norm( S_ref{n} ) / length(S_ref{n});
        Error_P(t,n-1) = norm( (P_ref{n} - P{n}) ) / length(P_ref{n});
        Error_S(t,n-1) = norm( (S_ref{n} - S{n}) ) / length(S_ref{n});
    end
        
%     figure(1); 
%     surf(xcm{2},ycm{2},reshape(abs(S_ref{2}-S{2}),Nx(2),Ny(2))');
%     view([0 0]);
%     title(strcat('|{S_w}^{DNS\_',num2str(Nx(1)),'\times',num2str(Ny(1)),'} - {S_w}^{DNS\_',num2str(Nx(2)),'\times',num2str(Ny(2)),'}|'));
%     xlabel('x[m]');
%     ylabel('y[m]');
%     hold on;
%     plot3([2:7],4.5*ones(6,1),0.0005*ones(6,1),'-r','LineWidth',3);
%     hold off;
end

% %% Ploting Pressure Error
% figure(2);
% legendInfo = cell(length(Nx)-1,1);
% for n=1:length(Nx)-1
%     legendInfo{n} = strcat('DNS\_',num2str(Nx(1)),'\times',num2str(Ny(1)),' vs EDFM\_',num2str(Nx(n+1)),'\times',num2str(Ny(n+1)));
%     plot( (1:t)*dt , Error_P(:,n) ,'LineWidth',2 );
%     hold on;
% end
% title('Pressure Error');
% legend(legendInfo);
% xlabel('Simulation Time [day], (\Deltat = 10^{-4} day)');
% ylabel('Error(P)');
% 
% %% Ploting Saturation Error
% figure(3);
% legendInfo = cell(length(Nx)-1,1);
% for n=1:length(Nx)-1
%     legendInfo{n} = strcat('DNS\_',num2str(Nx(1)),'\times',num2str(Ny(1)),' vs EDFM\_',num2str(Nx(n+1)),'\times',num2str(Ny(n+1)));
%     plot( (1:t)*dt , Error_S(:,n) ,'LineWidth',2 );
%     hold on;
% end
% title('Saturation Error');
% legend(legendInfo);
% xlabel('Simulation Time [day], (\Deltat = 10^{-4} day)');
% ylabel('Error(S)');
%% Manual Plot for 3 comparisons
figure(2);
hold on;
plot( (1:t)*dt , Error_P(:,1) , '-'  , 'LineWidth',2 );
plot( (1:t)*dt , Error_P(:,2) , '--' , 'LineWidth',2 );
plot( (1:t)*dt , Error_P(:,3) , '-.' , 'LineWidth',2 );
title('Pressure Error');
legend('DNS(225\times225) vs EDFM(45\times45)' , 'DNS(225\times225) vs EDFM(25\times25)' , 'DNS(225\times225) vs EDFM(15\times15)');
xlabel('Simulation Time [day], (\Deltat = 10^{-4} day)');
ylabel('Error(P)');

%% Manual Plot for 3 comparisons
figure(3);
hold on;
plot( (1:t)*dt , Error_S(:,1) , '-'  , 'LineWidth',2 );
plot( (1:t)*dt , Error_S(:,2) , '--' , 'LineWidth',2 );
plot( (1:t)*dt , Error_S(:,3) , '-.' , 'LineWidth',2 );
title('Saturation Error');
legend('DNS(225\times225) vs EDFM(45\times45)' , 'DNS(225\times225) vs EDFM(25\times25)' , 'DNS(225\times225) vs EDFM(15\times15)');
xlabel('Simulation Time [day], (\Deltat = 10^{-4} day)');
ylabel('Error(S)');
