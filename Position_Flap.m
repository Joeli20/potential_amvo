clear;
clc;
close all;

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

data_airfoil = "Airfoil_data_files/NACA_0010_N_"+N+"_coord.txt";
file_airfoil = fullfile(data_airfoil);

data_flap = "Airfoil_data_files/NACA_0015_N_"+N+"_coord.txt";
file_flap = fullfile(data_flap);

% Preallocating for speed
longitud_panel_airfoil = zeros(N,1);
control_points_airfoil = zeros(2,N);
cosinus_airfoil = zeros(1,N);
sinus_airfoil = zeros(1,N);
vector_normal_airfoil = zeros(2,N);
vector_tangent_airfoil = zeros(2,N);

longitud_panel_flap = zeros(N,1);
control_points_flap = zeros(2,N);
cosinus_flap = zeros(1,N);
sinus_flap = zeros(1,N);
vector_normal_flap = zeros(2,N);
vector_tangent_flap = zeros(2,N);

% Read the geometry coordinates from the text file
data_airfoil = importdata(file_airfoil);
nodal_points_airfoil(1,:) = data_airfoil(:,2);
nodal_points_airfoil(2,:) = data_airfoil(:,3);

data_flap = importdata(file_flap);
nodal_points_flap(1,:) = data_flap(:,2);
nodal_points_flap(2,:) = data_flap(:,3);

c_airfoil = 1;
c_flap = 0.45;
gap = 0.05; %Distance between trailing edge (airfoil) and leading edge (flap)
delta_flap = 45; %In degrees

% Escalar
nodal_points_airfoil(:,:) = nodal_points_airfoil(:,:)*c_airfoil;
nodal_points_flap_0(:,:) = nodal_points_flap(:,:)*c_flap;

% Posicionar + Rotar flap
nodal_points_flap(2,:) = nodal_points_flap_0(2,:)*cosd(delta_flap) - sind(delta_flap)*(gap + nodal_points_flap_0(1,:)); % Y component
nodal_points_flap(1,:) = c_airfoil + cosd(delta_flap)*(gap + nodal_points_flap_0(1,:)) + sind(delta_flap)*nodal_points_flap_0(2,:);% X component

% Calculations AIRFOIL

for i= 1:(length(nodal_points_airfoil)-1)
        longitud_panel_airfoil(i,1) = sqrt((nodal_points_airfoil(1,i+1)-nodal_points_airfoil(1,i))^2+(nodal_points_airfoil(2,i)-nodal_points_airfoil(2,i+1))^2);
        control_points_airfoil(1,i) = (nodal_points_airfoil(1,i+1)+nodal_points_airfoil(1,i))*0.5;
        control_points_airfoil(2,i) = (nodal_points_airfoil(2,i)+nodal_points_airfoil(2,i+1))*0.5;
        cosinus_airfoil(1,i) = (nodal_points_airfoil(1,i+1)-nodal_points_airfoil(1,i))/longitud_panel_airfoil(i,1);
        sinus_airfoil(1,i) = (nodal_points_airfoil(2,i)-nodal_points_airfoil(2,i+1))/longitud_panel_airfoil(i,1);
        vector_normal_airfoil(1,i) = sinus_airfoil(1,i);
        vector_normal_airfoil(2,i) = cosinus_airfoil(1,i);
        vector_tangent_airfoil(1,i) = cosinus_airfoil(1,i);
        vector_tangent_airfoil(2,i) = -sinus_airfoil(1,i);
end


% Plot the geometry of the surface with nodal and control points
figure;
plot(nodal_points_airfoil(1, [1:end, 1]), nodal_points_airfoil(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control_points_airfoil(1, :), control_points_airfoil(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Cylinder Geometry with Nodal and Control Points');
xlabel('X');
ylabel('Y');
legend('Nodal Points', 'Control Points');

% Optionally, plot the geometry along with the middle points
figure;
plot(nodal_points_airfoil(1, [1:end, 1]), nodal_points_airfoil(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control_points_airfoil(1, :), control_points_airfoil(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Geometria parametritzada');
xlabel('x');
ylabel('z');
legend('Node', 'Punt mig');
quiver(control_points_airfoil(1,:),control_points_airfoil(2,:),vector_normal_airfoil(1,:),vector_normal_airfoil(2,:));

% Calculations FLAP

for i= 1:(length(nodal_points_flap)-1)
        longitud_panel_flap(i,1) = sqrt((nodal_points_flap(1,i+1)-nodal_points_flap(1,i))^2+(nodal_points_flap(2,i)-nodal_points_flap(2,i+1))^2);
        control_points_flap(1,i) = (nodal_points_flap(1,i+1)+nodal_points_flap(1,i))*0.5;
        control_points_flap(2,i) = (nodal_points_flap(2,i)+nodal_points_flap(2,i+1))*0.5;
        cosinus_flap(1,i) = (nodal_points_flap(1,i+1)-nodal_points_flap(1,i))/longitud_panel_flap(i,1);
        sinus_flap(1,i) = (nodal_points_flap(2,i)-nodal_points_flap(2,i+1))/longitud_panel_flap(i,1);
        vector_normal_flap(1,i) = sinus_flap(1,i);
        vector_normal_flap(2,i) = cosinus_flap(1,i);
        vector_tangent_flap(1,i) = cosinus_flap(1,i);
        vector_tangent_flap(2,i) = -sinus_flap(1,i);
end


% Plot the geometry of the surface with nodal and control points
figure;
plot(nodal_points_airfoil(1, [1:end, 1]), nodal_points_airfoil(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control_points_airfoil(1, :), control_points_airfoil(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(nodal_points_flap(1, [1:end, 1]), nodal_points_flap(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
plot(control_points_flap(1, :), control_points_flap(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Airfoil + Flap Geometry with Nodal and Control Points');
xlabel('X');
ylabel('Y');
legend('Nodal Points', 'Control Points');

% Optionally, plot the geometry along with the middle points
figure;
plot(nodal_points_airfoil(1, [1:end, 1]), nodal_points_airfoil(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control_points_airfoil(1, :), control_points_airfoil(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(nodal_points_flap(1, [1:end, 1]), nodal_points_flap(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
plot(control_points_flap(1, :), control_points_flap(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Geometria parametritzada');
xlabel('x');
ylabel('z');
legend('Node', 'Punt mig');
quiver(control_points_airfoil(1,:),control_points_airfoil(2,:),vector_normal_airfoil(1,:),vector_normal_airfoil(2,:));
quiver(control_points_flap(1,:),control_points_flap(2,:),vector_normal_flap(1,:),vector_normal_flap(2,:));