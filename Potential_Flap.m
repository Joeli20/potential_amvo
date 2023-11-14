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

data_m = "Airfoil_data_files/NACA_0010_N_"+N+"_coord.txt";
file_m = fullfile(data_m); % Subindex _m means Main Airfoil

data_f = "Airfoil_data_files/NACA_0015_N_"+N+"_coord.txt";
file_f = fullfile(data_f); % Subindex _f means Flap Airfoil

% Read the geometry coordinates from the text file
data_m = importdata(file_m);
node_m(1,:) = data_m(:,2);
node_m(2,:) = data_m(:,3);

data_f = importdata(file_f);
node_f(1,:) = data_f(:,2);
node_f(2,:) = data_f(:,3);

clear data_m; clear data_f; % Maintaining Workspace clean
clear file_m; clear file_f; % Maintaining Workspace clean

%% PREALLOCATING
l_p_m = zeros(N,1);
control_m = zeros(2,N);
cosinus_m = zeros(1,N);
sinus_m = zeros(1,N);
vec_n_m = zeros(2,N);
vec_t_m = zeros(2,N);

l_p_f = zeros(N,1);
control_f = zeros(2,N);
cosinus_f = zeros(1,N);
sinus_f = zeros(1,N);
vec_n_f = zeros(2,N);
vec_t_f = zeros(2,N);

%% PARAMETERS DEFINITION
% Geometry
c_m = 1;
c_f = 0.45;
gap = 0.05; % distance between trailing edge (airfoil) and leading edge (flap)
delta_f = 45; % in degrees

% Fluid
Q_inf = 1; % in m/s
AoA = 0; % in degrees

%% GEOMETRY PREVIOUS CALCULATIONS
% Scalate
node_m(:,:) = node_m(:,:)*c_m;
node_f_0(:,:) = node_f(:,:)*c_f;

% Flap positioning + rotation
    % Y component
node_f(2,:) = node_f_0(2,:)*cosd(delta_f) - sind(delta_f)*(gap + node_f_0(1,:));
    % X component
node_f(1,:) = c_m + cosd(delta_f)*(gap + node_f_0(1,:)) + sind(delta_f)*node_f_0(2,:);

clear node_f_0;

%% AIRFOIL

for i= 1:(length(node_m)-1)
        l_p_m(i,1) = sqrt((node_m(1,i+1)-node_m(1,i))^2+(node_m(2,i)-node_m(2,i+1))^2);
        control_m(1,i) = (node_m(1,i+1)+node_m(1,i))*0.5;
        control_m(2,i) = (node_m(2,i)+node_m(2,i+1))*0.5;
        cosinus_m(1,i) = (node_m(1,i+1)-node_m(1,i))/l_p_m(i,1);
        sinus_m(1,i) = (node_m(2,i)-node_m(2,i+1))/l_p_m(i,1);
        vec_n_m(1,i) = sinus_m(1,i);
        vec_n_m(2,i) = cosinus_m(1,i);
        vec_t_m(1,i) = cosinus_m(1,i);
        vec_t_m(2,i) = -sinus_m(1,i);
end
clear i;

%% FLAP

for i= 1:(length(node_f)-1)
        l_p_f(i,1) = sqrt((node_f(1,i+1)-node_f(1,i))^2+(node_f(2,i)-node_f(2,i+1))^2);
        control_f(1,i) = (node_f(1,i+1)+node_f(1,i))*0.5;
        control_f(2,i) = (node_f(2,i)+node_f(2,i+1))*0.5;
        cosinus_f(1,i) = (node_f(1,i+1)-node_f(1,i))/l_p_f(i,1);
        sinus_f(1,i) = (node_f(2,i)-node_f(2,i+1))/l_p_f(i,1);
        vec_n_f(1,i) = sinus_f(1,i);
        vec_n_f(2,i) = cosinus_f(1,i);
        vec_t_f(1,i) = cosinus_f(1,i);
        vec_t_f(2,i) = -sinus_f(1,i);
end
clear i;

%% PLOTTING

% Airfoil geometry
figure;
plot(node_m(1, [1:end, 1]), node_m(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control_m(1, :), control_m(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Cylinder Geometry with Nodal and Control Points');
xlabel('X');
ylabel('Y');
legend('Nodal Points', 'Control Points');

% Airfoil vectors
figure;
plot(node_m(1, [1:end, 1]), node_m(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control_m(1, :), control_m(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Geometria parametritzada');
xlabel('x');
ylabel('z');
legend('Node', 'Punt mig');
quiver(control_m(1,:),control_m(2,:),vec_n_m(1,:),vec_n_m(2,:));

% Airfoil + Flap geometry
figure;
plot(node_m(1, [1:end, 1]), node_m(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control_m(1, :), control_m(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(node_f(1, [1:end, 1]), node_f(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
plot(control_f(1, :), control_f(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Airfoil + Flap Geometry with Nodal and Control Points');
xlabel('X');
ylabel('Y');
legend('Nodal Points', 'Control Points');

% Airfoil + Flap vectors
figure;
plot(node_m(1, [1:end, 1]), node_m(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control_m(1, :), control_m(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(node_f(1, [1:end, 1]), node_f(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
plot(control_f(1, :), control_f(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Geometria parametritzada');
xlabel('x');
ylabel('z');
legend('Node', 'Punt mig');
quiver(control_m(1,:),control_m(2,:),vec_n_m(1,:),vec_n_m(2,:));
quiver(control_f(1,:),control_f(2,:),vec_n_f(1,:),vec_n_f(2,:));