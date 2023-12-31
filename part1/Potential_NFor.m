clear;
clc;
close all;

%% DATA INPUT
% Specify the file path for the geometry coordinates text file
for k = 1:6
clear node; clear data;

N = 2^(k+3);

data = "Airfoil_data_files/NACA_0010_N_"+N+"_coord.txt";
file = fullfile(data); % Subindex _m means Main Airfoil

% Read the geometry coordinates from the text file
data = importdata(file);
node(1,:) = data(:,2);
node(2,:) = data(:,3);

clear data; clear file; % Maintaining Workspace clean

%% PREALLOCATING
l_p = zeros(N,1);
control = zeros(2,N);
cosinus = zeros(1,N);
sinus = zeros(1,N);
vec_n = zeros(2,N);
vec_t = zeros(2,N);

%% PARAMETERS DEFINITION
% Geometry
c = 1;

%% GEOMETRY PREVIOUS CALCULATIONS
% Scalate
node(:,:) = node(:,:)*c;

%% AIRFOIL

for i= 1:(length(node)-1)
        l_p(i,1) = sqrt((node(1,i+1)-node(1,i))^2+(node(2,i)-node(2,i+1))^2);
        control(1,i) = (node(1,i+1)+node(1,i))*0.5;
        control(2,i) = (node(2,i)+node(2,i+1))*0.5;
        cosinus(1,i) = (node(1,i+1)-node(1,i))/l_p(i,1);
        sinus(1,i) = (node(2,i)-node(2,i+1))/l_p(i,1);
        vec_n(1,i) = sinus(1,i);
        vec_n(2,i) = cosinus(1,i);
        vec_t(1,i) = cosinus(1,i);
        vec_t(2,i) = -sinus(1,i);
end
clear i;

%% APARTAT 1
Q_inf = 1; % in m/s
AoA_1 = [2 4 6 8 10]; % in degrees

% Preallocating
v_f_1 = zeros(length(AoA_1),N);
v_x_1 = zeros(length(AoA_1),N);
v_z_1 = zeros(length(AoA_1),N);
cp_1 = zeros(length(AoA_1),N);
gamma_1 = zeros(length(AoA_1),N);

for i = 1:length(AoA_1)
    [v_f,v_x,v_z,cp,cl,cm_0,gamma] = Vortex(Q_inf,AoA_1(i),cosinus,sinus,l_p,node,control,vec_t,c);
    v_f_1(i,:) = v_f;
    v_x_1(i,:) = v_x;
    v_z_1(i,:) = v_z;
    cp_1(i,:) = cp;
    cl_1(i,k) = cl;
    cm_0_1(i,k) = cm_0;
    gamma_1(i,:) = gamma;
end
clear i;

%% APARTAT 2
AoA_2 = [0 2 4 6];
M_inf = linspace(0.01, 0.99, 100);
GAMMA = 1.4;

% Preallocating
v_f_2 = zeros(N,length(M_inf),length(AoA_2));
v_x_2 = zeros(N,length(M_inf),length(AoA_2));
v_z_2 = zeros(N,length(M_inf),length(AoA_2));
cp_2 = zeros(N,length(M_inf),length(AoA_2));
cl_2 = zeros(length(M_inf),length(AoA_2));
cm_0_2 = zeros(length(M_inf),length(AoA_2));
gamma_2 = zeros(N,length(M_inf),length(AoA_2));

for i = 1:length(AoA_2)
    for j = 1:length(M_inf)
        [v_f,v_x,v_z,cp,cl,cm_0,gamma] = Vortex(Q_inf,AoA_2(i),cosinus,sinus,l_p,node,control,vec_t,c);
        v_f_2(:,j,i) = v_f;
        v_x_2(:,j,i) = v_x;
        v_z_2(:,j,i) = v_z;
        cp_2(:,j,i) = cp;
        cl_2(j,i) = cl;
        cm_0_2(j,i) = cm_0;
        gamma_2(:,j,i) = gamma;

        cp_min(i) = min(cp_2(:,j,i));

        cp_star(j,i) = (2/(GAMMA*M_inf(j).^2))*(((2+(GAMMA-1)*M_inf(j).^2)/(1+GAMMA))^(GAMMA/(GAMMA-1))-1);
        cp_Laitone(j,i) = cp_min(i)/(sqrt(1-M_inf(j)^2)+(cp_min(i)/2)*(M_inf(j)^2/sqrt(1-M_inf(j)^2))*(1+(GAMMA-1)*M_inf(j)^2/2));

        if abs(cp_Laitone(j,i)-cp_star(j,i))<0.1
            M_inf_2(i,k) = M_inf(j);
            break;
        end
    end
end
clear cp_min; clear cp_star; clear cp_Laitone; clear i;

%% APARTAT 3
AoA_3 = 2;
M_inf_3 = [M_inf_2(2)-0.15 M_inf_2(2)-0.1 M_inf_2(2)-0.05 M_inf_2(2)];

% Preallocating
v_f_3 = zeros(length(M_inf_3),N);
v_x_3 = zeros(length(M_inf_3),N);
v_z_3 = zeros(length(M_inf_3),N);
cp_3 = zeros(length(M_inf_3),N);
cl_3_incompressible = zeros(1,length(M_inf_3));
cm_0_3 = zeros(1,length(M_inf_3));
gamma_3 = zeros(length(M_inf_3),N);


for i = 1:length(M_inf_3)
    Beta(i) = sqrt(1-M_inf_3(i)^2);
    [v_f,v_x,v_z,cp,cl,cm_0,gamma] = Vortex(Q_inf,AoA_3,cosinus,sinus,l_p,node,control,vec_t,c);
    v_f_3(i,:) = v_f;
    v_x_3(i,:) = v_x;
    v_z_3(i,:) = v_z;
    cp_3(i,:) = cp;
    cl_3_incompressible(i) = cl;
    cm_0_3(i) = cm_0;
    gamma_3(i,:) = gamma;
    cl_3(i,k) = cl_3_incompressible(i)./Beta(i);
end
end

%% PLOTTING GENERAL

% Airfoil geometry
figure;
plot(node(1, [1:end, 1]), node(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control(1, :), control(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Cylinder Geometry with Nodal and Control Points');
xlabel('X');
ylabel('Y');
legend('Nodal Points', 'Control Points');

% Airfoil vectors
figure;
plot(node(1, [1:end, 1]), node(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control(1, :), control(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Geometria parametritzada');
xlabel('x');
ylabel('z');
legend('Node', 'Punt mig');
quiver(control(1,:),control(2,:),vec_n(1,:),vec_n(2,:));

%% PLOTTING PART 1
% Cl vs AoA
figure;
plot(AoA_1,cl_1, '-', 'LineWidth', 2);
title('Cl vs AoA for NACA 0010');
xlabel('Angle of attack');
ylabel('Cl');
legend('Cl');

% Cm vs AoA
figure;
plot(AoA_1,cm_0_1, '-', 'LineWidth', 2);
title('Cm_1_/_4 vs AoA for NACA 0010');
xlabel('Angle of attack');
ylabel('Cm_1_/_4');
legend('Cm_1_/_4');

%% PLOTTING PART 2
% Cl vs AoA
AoA_2 = [0 2 4 6];
figure;
plot(AoA_2,M_inf_2, '-', 'LineWidth', 2);
title('M_c_r vs AoA for NACA 0010');
xlabel('Angle of attack');
ylabel('M_c_r');
legend('M_c_r');