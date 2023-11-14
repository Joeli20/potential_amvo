clear;
clc;
close all;

%% DATA INPUT
% Specify the file path for the geometry coordinates text file
e = 0;
while e~=1
panells = "Number of divisions (16,32,64,128,256,512):";
N = input(panells);

if N~=16 && N~=32 && N~=64 && N~=128 && N~=256 && N~=512
    clc;
    disp("Number of divisions not valid. Please enter a correct value.");
else
    e = 1;
end
end

clear e; clear panells; % Maintaining Workspace clean

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

%[v_f,cp,a_ii,sigma] = Sources(Q_inf,AoA,cosinus,sinus,l_p,node,control,vec_n);
%% APARTAT 1
Q_inf = 1; % in m/s
AoA_1 = [2 4 6 8 10]; % in degrees

% Preallocating
v_f_1 = zeros(length(AoA_1),N);
v_x_1 = zeros(length(AoA_1),N);
v_z_1 = zeros(length(AoA_1),N);
cp_1 = zeros(length(AoA_1),N);
cl_1 = zeros(1,length(AoA_1));
cm_0_1 = zeros(1,length(AoA_1));
gamma_1 = zeros(length(AoA_1),N);

for i = 1:length(AoA_1)
[v_f,v_x,v_z,cp,cl,cm_0,gamma] = Vortex(Q_inf,AoA_1(i),cosinus,sinus,l_p,node,control,vec_t,c);
v_f_1(i,:) = v_f;
v_x_1(i,:) = v_x;
v_z_1(i,:) = v_z;
cp_1(i,:) = cp;
cl_1(i) = cl;
cm_0_1(i) = cm_0;
gamma_1(i,:) = gamma;
end
clear i;

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