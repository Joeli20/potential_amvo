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

clear e; clear panells; % Maintaining Workspace Clean

data = "Airfoil_data_files/NACA_0010_N_"+N+"_coord.txt";
file = fullfile(data);

% Preallocating for speed
l_p = zeros(N,1);
control = zeros(2,N);
cosinus = zeros(1,N);
sinus = zeros(1,N);
vec_n = zeros(2,N);
vec_t = zeros(2,N);

% Read the geometry coordinates from the text file
data = importdata(file);

node(1,:) = data(:,2);
node(2,:) = data(:,3);

clear data; clear file;

% Initial conditions
Q_inf = 1;
AoA = 0;

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

[v_f,cp,a_ii,sigma] = Sources(Q_inf,cosinus,sinus,l_p,node,control,vec_n);

% Plot the geometry of the cylinder with nodal and control points
figure;
plot(node(1, [1:end, 1]), node(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control(1, :), control(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Cylinder Geometry with Nodal and Control Points');
xlabel('X');
ylabel('Y');
legend('Nodal Points', 'Control Points');

% Optionally, plot the geometry along with the middle points
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