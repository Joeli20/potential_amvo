clear all; clc;

% Specify the file path for the geometry coordinates text file
file = 'naca_parameters.txt';

% Read the geometry coordinates from the text file
data = importdata(file);

nodal_points(1,:) = data(:,1);
nodal_points(2,:) = -data(:,2);

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