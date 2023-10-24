% Define cylinder parameters
cylinder_radius = 1.0;
cylinder_center = [0, 0];
num_panels = 10; % Number of panels

% Create parametric values for the cylinder
theta = linspace(0, 2 * pi, num_panels + 1); % Angles

% Nodal and control point arrays (including an extra point to close the circle)
nodal_points2 = [cylinder_radius * cos(theta); cylinder_radius * sin(theta)];
control_points2 = 0.5 * (nodal_points2(:, 1:end-1) + nodal_points2(:, 2:end));

% Calculate the panel lengths
panel_lengths = sqrt(diff(nodal_points(1, :)).^2 + diff(nodal_points(2, :)).^2);

% Calculate normal vectors
normal_vectors = [-diff(nodal_points(2, :)); diff(nodal_points(1, :))];
normal_vectors = normal_vectors ./ vecnorm(normal_vectors);

% Calculate tangent vectors
tangent_vectors = [diff(nodal_points(1, :)); diff(nodal_points(2, :))];
tangent_vectors = tangent_vectors ./ vecnorm(tangent_vectors);


% Plot the geometry of the cylinder with nodal and control points
figure;
plot(nodal_points2(1, [1:end, 1]), nodal_points2(2, [1:end, 1]), 'bo-', 'LineWidth', 2);
hold on;
plot(control_points2(1, :), control_points2(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Cylinder Geometry with Nodal and Control Points');
xlabel('X');
ylabel('Y');
legend('Nodal Points', 'Control Points');
