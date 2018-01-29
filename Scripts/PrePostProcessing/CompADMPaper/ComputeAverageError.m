% Folders Molar ADM
clc
Directory_ADM = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_ADM_natural/';
Directory_FS = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_FS_natural/';
Angle = ['00deg'; '15deg'; '30deg'; '45deg'; '98deg'];
dz = ['0.05';'0.07';'0.10'];
Error_P_norm = zeros(5, 3, 20);
Error_S_norm = zeros(5, 3, 20);
Av_nodes = zeros(5, 3, 20);
max_nodes = zeros(5, 3, 20);
min_nodes = zeros(5, 3, 20);
MaxTimeStep=10;

for i = 1:5
    disp(Angle(i,:));
    for j=1:20
        %% Directory FS Solutions
        Dir_fs = strcat(Directory_FS,Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/Output/Solution/');
        Name_fs = strcat('Case1_',Angle(i,:),'_', num2str(j),'_FS_natural_Sol');
        for k=1:3
            %% Directory ADM Solutions
            Dir_adm = strcat(Directory_ADM,Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/dz_',num2str(k),'/Output/Solution/');
            Name_adm = strcat('Case1_',Angle(i,:),'_', num2str(j),'natural_ADM_dz_', num2str(k),'_Sol');
            Dir_nodes = strcat(Directory_ADM,Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/dz_',num2str(k),'/Output/');
            Name_nodes = 'ADMStats.txt';
            Nodes = load(strcat(Dir_nodes, Name_nodes));
            Error_P = zeros(10, 1);
            Error_S = zeros(10, 1);
            for time = 1:MaxTimeStep
                % Error norm for time t
                Solution_FS = load(strcat(Dir_fs, Name_fs, num2str(time),'.txt'));
                Solution_adm = load(strcat(Dir_adm, Name_adm, num2str(time),'.txt'));
                Error_P(time) = norm(Solution_FS(:,2) - Solution_adm(:,2), 2)/100;
                Error_S(time) = mean(abs(Solution_FS(:,3) - Solution_adm(:,3)));
            end
            % Average over times
            Error_P_norm(i, k, j) = mean(Error_P);
            Error_S_norm(i, k, j) = mean(Error_S);
            
            % AVerage nodes for dz_i
            Av_nodes(i, k, j) = mean(Nodes(:,4));
            max_nodes(i, k, j) = max(Nodes(:, 4));
            min_nodes(i, k, j) = min(Nodes(:, 4));
        end
    end
end
disp(char(5))

%% Average error
for i=1:5
    disp(Angle(i,:));
    for k=1:3
        disp(['dz_', num2str(k)])
        AV_error_angle = mean(Error_P_norm(i,k,:));
        disp('Pressure')
        disp(AV_error_angle);
        AV_error_angle = mean(Error_S_norm(i,k,:));
        disp('Saturation')
        disp(AV_error_angle);
    end
end

for i=1:5
    disp(Angle(i,:));
     for k=1:3
        disp(['dz_', num2str(k)])
        av = mean(Av_nodes(i,k,:));
        max = mean(max_nodes(i,k,:));
        min = mean(min_nodes(i,k,:));
        disp([av, min, max]);
     end
end