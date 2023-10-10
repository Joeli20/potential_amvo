% Specify the file path for the geometry coordinates text file
file_path = 'naca_parameters.txt';

% Read the geometry coordinates from the text file
geometry_data = importdata(file_path);

% Extract X and Y coordinates
X = geometry_data(:, 1);
Y = geometry_data(:, 2);

% Calculate the number of points
num_points = length(X);

% Connect the first and last points to close the geometry
X = [X; X(1)];
Y = [Y; Y(1)];

% Calculate the coordinates of nodal and control points
nodal_points = [X'; Y']; % Nodal points are now closed
control_points = 0.5 * (nodal_points(:, 1:end-1) + nodal_points(:, 2:end));

% Calculate the middle points of each panel
mid_points = 0.5 * (nodal_points(:, 1:end-1) + nodal_points(:, 2:end));

% Calculate the panel lengths
panel_lengths = sqrt(diff(nodal_points(1, :)).^2 + diff(nodal_points(2, :)).^2);

% Calculate the number of subdivisions (panels)
num_subdivisions = length(panel_lengths);

% Calculate normal vectors
normal_vectors = [-diff(nodal_points(2, :)); diff(nodal_points(1, :))];
normal_vectors = normal_vectors ./ vecnorm(normal_vectors);

% Calculate tangent vectors
tangent_vectors = [diff(nodal_points(1, :)); diff(nodal_points(2, :))];
tangent_vectors = tangent_vectors ./ vecnorm(tangent_vectors);

% Optionally, plot the geometry along with the middle points
figure;
plot(X, Y, 'bo-', 'LineWidth', 2);
hold on;
plot(mid_points(1, :), mid_points(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Parametric Geometry with Middle Points');
xlabel('X');
ylabel('Y');
legend('Geometry Points', 'Middle Points');
