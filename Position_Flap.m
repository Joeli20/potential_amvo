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
longitud_panel = zeros(N,1);
control_points = zeros(2,N);
cosinus = zeros(1,N);
sinus = zeros(1,N);
vector_normal = zeros(2,N);
vector_tangent = zeros(2,N);
nodal_points = zeros(2,(N+1)*2);

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
delta_flap = 32; %In degrees

% Escalar
nodal_points_airfoil(:,:) = nodal_points_airfoil(:,:)*c_airfoil;
nodal_points_flap(:,:) = nodal_points_flap(:,:)*c_flap;

% Posicionar + Rotar flap
nodal_points_flap(2,:) = nodal_points_flap(2,:) - sind(delta_flap)*(gap + nodal_points_flap(1,:)); % Y component
nodal_points_flap(1,:) = c_airfoil + cosd(delta_flap)*(gap + nodal_points_flap(1,:));% X component

nodal_points = [nodal_points_airfoil,nodal_points_flap];

% Calculations

for i= 1:(length(nodal_points)-1)
        longitud_panel(i,1) = sqrt((nodal_points(1,i+1)-nodal_points(1,i))^2+(nodal_points(2,i)-nodal_points(2,i+1))^2);
        control_points(1,i) = (nodal_points(1,i+1)+nodal_points(1,i))*0.5;
        control_points(2,i) = (nodal_points(2,i)+nodal_points(2,i+1))*0.5;
        cosinus(1,i) = (nodal_points(1,i+1)-nodal_points(1,i))/longitud_panel(i,1);
        sinus(1,i) = (nodal_points(2,i)-nodal_points(2,i+1))/longitud_panel(i,1);
        vector_normal(1,i) = sinus(1,i);
        vector_normal(2,i) = cosinus(1,i);
        vector_tangent(1,i) = cosinus(1,i);
        vector_tangent(2,i) = -sinus(1,i);
end


% Plot the geometry of the cylinder with nodal and control points
figure;
plot(nodal_points(1, [1:end, 1]), nodal_points(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control_points(1, :), control_points(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Cylinder Geometry with Nodal and Control Points');
xlabel('X');
ylabel('Y');
legend('Nodal Points', 'Control Points');

% Optionally, plot the geometry along with the middle points
figure;
plot(nodal_points(1, [1:end, 1]), nodal_points(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control_points(1, :), control_points(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Geometria parametritzada');
xlabel('x');
ylabel('z');
legend('Node', 'Punt mig');
quiver(control_points(1,:),control_points(2,:),vector_normal(1,:),vector_normal(2,:));



% Integral of convective term, x and y components.
%
% Written by: Joel Campo, Jordi Gallart, Martí Santamaria, 2023
% Group 16. AMVO. MUEA.
%
% Inputs:
%   u: Matrix of the horitzontal velocity components
%   v: Matrix of the vertical velocity components
%   L: length of a side of the analysed square
% Outputs:
%   u_conv_num: Solution of the convective terms of horitzontal velocity
%   v_conv_num: Solution of the convective terms of vertical velocity