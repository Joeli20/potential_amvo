clear; clc; close all;

% Define cylinder parameters
cylinder_radius = 1.0;
cylinder_center = [0, 0];
N = 32; % Number of panels

% Create parametric values for the cylinder
theta = linspace(0, 2 * pi, N + 1); % Angles

% Nodal arrays
nodal_points = [cylinder_radius * cos(theta); -cylinder_radius * sin(theta)];

% Preallocating for speed
longitud_panel = zeros(N,1);
control_points = zeros(2,N);
cosinus = zeros(1,N);
sinus = zeros(1,N);
vector_normal = zeros(2,N);
vector_tangent = zeros(2,N);

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

[V_final,Cp, a_ii, sigma] = ViP_fons(1,cosinus,sinus,longitud_panel,nodal_points,control_points,vector_normal);
[V_final_vortex,Cp_vortex, a_ii_vortex, gamma] = ViP_vortex(1,cosinus,sinus,longitud_panel,nodal_points,control_points,vector_normal, vector_tangent);
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

angle_distribution = zeros(1,N);
for i = 1:N
angle_distribution(1,i)= (360/N)*(i-1);
end

figure;
plot(angle_distribution,sigma);
title('Sigma vs angle');
xlabel('angle º');
ylabel('sigma');
legend('sigma');

figure;
plot(angle_distribution,V_final);
title('Velocity vs angle');
xlabel('angle º');
ylabel('Velocity');
legend('Velocity');

figure;
plot(angle_distribution,Cp);
title('Cp vs angle');
xlabel('angle º');
ylabel('Cp');
legend('Cp');