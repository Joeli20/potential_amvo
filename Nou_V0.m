
clear all; clc;

% Define cylinder parameters
cylinder_radius = 1.0;
cylinder_center = [0, 0];
num_panels = 20; % Number of panels

% Create parametric values for the cylinder
theta = linspace(0, 2 * pi, num_panels + 1); % Angles

% Nodal arrays
nodal_points = [cylinder_radius * cos(theta); -cylinder_radius * sin(theta)];

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

[V_final,Cp, a_ii, sigma, V] = ViP_fons(10,cosinus,sinus,longitud_panel,nodal_points,control_points,vector_normal);


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