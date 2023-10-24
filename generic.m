% Specify the file path for the geometry coordinates text file
file = 'naca_parameters.txt';

% Read the geometry coordinates from the text file
data = importdata(file);

% Extract X and Y coordinates
x = data(:, 1);
z = data(:, 2);

% Calculate the number of points
num_points = length(x);

% Connect the first and last points to close the geometry
x = [x; x(1)];
z = [z; z(1)];

% Calculate the coordinates of nodal and control points
nodes = [x'; z']; % Nodal points are now closed

% Calculate the middle points of each panel
puntmig = 0.5 * (nodes(:, 1:end-1) + nodes(:, 2:end));

% Calculate the panel lengths
l_p = sqrt(diff(nodes(1, :)).^2 + diff(nodes(2, :)).^2);

% Calculate the number of subdivisions (panels)
N = length(l_p);

% Calculate normal vectors
normal = [-diff(nodes(2, :)); diff(nodes(1, :))];
normal = -normal ./ vecnorm(normal);

% Calculate tangent vectors
tangent = [diff(nodes(1, :)); diff(nodes(2, :))];
tangent = tangent ./ vecnorm(tangent);

for i = 1:N
    alpha_i(i) = atan2(normal(1,i),normal(2,i)); %Repassar perquè al refer conversio canvien signes
end

% Optionally, plot the geometry along with the middle points
figure;
plot(x, z, 'bo-', 'LineWidth', 2);
hold on;
plot(puntmig(1, :), puntmig(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Geometria parametritzada');
xlabel('x');
ylabel('z');
legend('Node', 'Punt mig');
quiver(puntmig(1,:),puntmig(2,:),normal(1,:),normal(2,:));
