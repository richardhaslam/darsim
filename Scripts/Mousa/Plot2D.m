%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Plot2D  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, TU Delft
% Project: F-ADM, 2016-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2017-12-04
% Modified on: 2017-12-04
%%
function Plot2D(DiscretizationModel,ProductionSystem,varType)
transparency = 'Transparent';
if DiscretizationModel.ReservoirGrid.Nz~=1
    error('This test case is not in 2D. Therefore "Plot2D" cannot provide visulization!\n');
end
NumOfFrac = ProductionSystem.FracturesNetwork.NumOfFrac;

% Reservoir Data
Nm = [DiscretizationModel.ReservoirGrid.Nx
      DiscretizationModel.ReservoirGrid.Ny];

Lm = [ProductionSystem.Reservoir.Length
      ProductionSystem.Reservoir.Width];

dm = Lm./Nm;
Xm = linspace( dm(1) , Lm(1)-dm(1) , Nm(1) )';
Ym = linspace( dm(2) , Lm(2)-dm(2) , Nm(2) )';

% Fractures Data
Nf = zeros(NumOfFrac,1);
Lf = zeros(NumOfFrac,1);
df = zeros(NumOfFrac,1);
Point_Start = zeros(2,NumOfFrac);
Point_End = zeros(2,NumOfFrac);
Xf = cell(NumOfFrac,1);
Yf = cell(NumOfFrac,1);
for f = 1 : NumOfFrac
    if DiscretizationModel.FracturesGrid.Grids(f).Ny~=1
        error('This test case has a fracture in 2D. Fractures should only be 1D lines. Therefore "Plot2D" cannot provide visulization!\n');
    end
    Nf(f) = DiscretizationModel.FracturesGrid.Grids(f).Nx;
    Lf(f) = ProductionSystem.FracturesNetwork.Fractures(f).Length;
    df(f) = Lf(f)/Nf(f);
    
    Point_Start(:,f) = [DiscretizationModel.FracturesGrid.Grids(f).GridCoords(1,1)
                        DiscretizationModel.FracturesGrid.Grids(f).GridCoords(1,2)];
    
    Point_End(:,f)   = [DiscretizationModel.FracturesGrid.Grids(f).GridCoords(end,1)
                        DiscretizationModel.FracturesGrid.Grids(f).GridCoords(end,2)];
    
    dfx = (Point_End(1,f)-Point_Start(1,f))/Nf(f);
    dfy = (Point_End(2,f)-Point_Start(2,f))/Nf(f);
    Xf{f} = linspace( Point_Start(1,f)+dfx , Point_End(1,f)-dfx , Nf(f) )';
    Yf{f} = linspace( Point_Start(2,f)+dfx , Point_End(2,f)-dfy , Nf(f) )';
end

Var = zeros( prod(Nm)+sum(Nf) , 1 );
Var(1:prod(Nm)) = ProductionSystem.Reservoir.State.Properties(varType).Value;
for f = 1 : NumOfFrac
    If_Start = prod(Nm)+sum(Nf(1:f-1))+1;
    If_End   = If_Start + Nf(f)-1;
    Var(If_Start:If_End) = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(varType).Value;
end

%% Plotting the variable
figure;
Z = reshape(Var(1:prod(Nm)),Nm(1),Nm(2))';
surf(Xm,Ym,Z);

%% Plotting Wells
Inj = ProductionSystem.Wells.Inj;
Prod = ProductionSystem.Wells.Prod;
hold on;
for w = 1 : length(Inj)
    X = Xm(Inj(w).Coord(1,1));
    Y = Ym(Inj(w).Coord(2,1));
    plot3(X,Y,Inj(w).p,'or','LineWidth',3);
end
for w = 1 : length(Prod)
    X = Xm(Prod(w).Coord(1,1));
    Y = Ym(Prod(w).Coord(2,1));
    plot3(X,Y,Prod(w).p,'ob','LineWidth',3);
end
hold off;

%% Plotting Fractures
hold on;
for f = 1 : NumOfFrac
    If_Start = prod(Nm)+sum(Nf(1:f-1))+1;
    If_End   = If_Start + Nf(f)-1;
    plot3(Xf{f},Yf{f},Var(If_Start:If_End)*1.01,'r','LineWidth',4);
%     if (length(frac) <= 5) && ( strcmp(func_name,'Error'  )==0 )
%         plot3(frac(f).xcf,frac(f).ycf,max(func)*ones(Nf(f),1),'r','LineWidth',4);
%     end
end
hold off;

%% Changing(reducing) the number of gridlines (edges / mesh) 
shading flat;
if     (Nm(1) <= 25)                   ,  alpha_v = 1.0;
elseif (Nm(1) >  25 ) && (Nm(1) <= 50 ),  alpha_v = 0.8;
elseif (Nm(1) >  50 ) && (Nm(1) <= 75 ),  alpha_v = 0.5;
elseif (Nm(1) >  75 ) && (Nm(1) <= 100),  alpha_v = 0.3;
elseif (Nm(1) > 100 )                  ,  alpha_v = 0.1;
end
if     (Nm(2) <= 25 )                  ,  alpha_h = 1.0;
elseif (Nm(2) >  25 ) && (Nm(2) <= 50 ),  alpha_h = 0.8;
elseif (Nm(2) >  50 ) && (Nm(2) <= 75 ),  alpha_h = 0.5;
elseif (Nm(2) >  75 ) && (Nm(2) <= 100),  alpha_h = 0.3;
elseif (Nm(2) > 100 )                  ,  alpha_h = 0.1;
end

Z_node = zeros( Nm(1)+1 , Nm(2)+1 );
for i = 1 : Nm(1)+1
    for j = 1 : Nm(2)+1
        if i>1
            if j>1      ,  Z1 = Z( i-1 , j-1 );  else,  Z1 = Z( i-1 , j   );  end
            if j<Nm(2)+1,  Z2 = Z( i-1 , j   );  else,  Z2 = Z( i-1 , j-1 );  end
        else
            if j>1      ,  Z1 = Z( i   , j-1 );  else,  Z1 = Z( i   , j   );  end
            if j<Nm(2)+1,  Z2 = Z( i   , j   );  else,  Z2 = Z( i   , j-1 );  end
        end
        if i<Nm(1)+1
            if j>1      ,  Z3 = Z( i   , j-1 );  else,  Z3 = Z( i   , j   );  end
            if j<Nm(2)+1,  Z4 = Z( i   , j   );  else,  Z4 = Z( i   , j-1 );  end
        else
            if j>1      ,  Z3 = Z( i-1 , j-1 );  else,  Z3 = Z( i-1 , j   );  end
            if j<Nm(2)+1,  Z4 = Z( i-1 , j   );  else,  Z4 = Z( i-1 , j-1 );  end
        end
        Z_node( i , j ) = (Z1+Z2+Z3+Z4)/4;
        Z_node( i , j ) = max([Z1 Z2 Z3 Z4])*1.001;
    end
end
hold on;
% Plotting lines in the X-Z plane
Xh = linspace( 0 , Lm(1) , Nm(1)+1 );
Yv = linspace( 0 , Lm(2) , Nm(2)+1 );
for j = 1 : Nm(2)+1
    Yh = ones ( Nm(2)+1 , 1) * (j-1) * dm(2);
    plot3(Xh,Yh,Z_node(j,:),'Color',[0 0 0 alpha_h]);
end
% Plotting lines in the Y-Z plane
for i = 1 : Nm(1)+1
    Xv = ones ( Nm(1)+1 , 1) * (i-1) * dm(1);
    plot3(Xv,Yv,Z_node(:,i),'Color',[0 0 0 alpha_v]);
end
hold off

%% Drawing gridlines for primal coarse cell boundaries
MSL = 0;
if isprop(DiscretizationModel,'Nc')
    MSL = DiscretizationModel.maxLevel(1);
    MSR = [DiscretizationModel.Coarsening(1,1,1);
           DiscretizationModel.Coarsening(1,2,1)];
    hold on
    for CL = 1:MSL
        ms_color = [1 1 1 0.3] - (CL-1)*0.05;
        if CL == MSL,  ms_color = [1 1 1 1] - (CL-1)*0.05;  end
        % Plotting lines in the X-Z plane
        for j = MSR(1)^CL+1:MSR(1)^CL:Nm(2)-MSR(1)^CL+1
            Yh = ones ( Nm(2)+1 , 1) * (j-1) * dm(2);
            plot3(Xh,Yh,Z_node(j,:)*(1+0.001*CL),'Color',ms_color,'LineWidth',2*CL-1);
        end
        % Plotting lines in the Y-Z plane
        for i = MSR(2)^CL+1:MSR(2)^CL:Nm(1)-MSR(2)^CL+1
            Xv = ones ( Nm(1)+1 , 1) * (i-1) * dm(1);
            plot3(Xv,Yv,Z_node(:,i)*(1+0.001*CL),'Color',ms_color,'LineWidth',2*CL-1);
        end
    end
    hold off   
end

%% Title, Axis, Colorbar and Rotation
switch varType
    case 'P_1'
        varText = 'Pressure';
    case 'S_1'
        varText = 'Saturation';
end
if MSL == 0
    scaleText = 'Finescale';
else
    scaleText = strcat(num2str(MSL), ' Level Multiscale');
end
title([scaleText, ' ', varText]);

set(0, 'DefaultAxesFontSize', 20, 'DefaultAxesFontWeight', 'normal');
axis([0 Lm(1) 0 Lm(2) min([Prod(:).p]) max([Inj(:).p]) ]);
if ( strcmp(transparency,'Transparent')==1 ),  alpha(0.8);  end
xlabel('x[m]');
ylabel('y[m]');

colorbar;
view(0, 90);
rotate3d on;
end