%%
clear; clc; close all;

% Define cylinder parameters
cylinder_radius = 1;
cylinder_center = [0, 0];
N = 199; % Number of panels

Qinf=2;
AoA=6;
c=2;

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
clear i;
%[V_final,Cp, a_ii, sigma] = Sources(1,cosinus,sinus,longitud_panel,nodal_points,control_points,vector_normal);
%[V_final_vortex, Vx, Vz, Cp_vortex,CL_perfil, C_l_cada_punt, Lift, CMO, Cm_0, a_ii_vortex, gamma] = Vortex(Qinf,AoA,cosinus,sinus,longitud_panel,nodal_points,control_points,vector_normal, vector_tangent);

%% APARTAT 2
AoA_2 = [0 2 4 6];
M_inf = linspace(0.01, 0.99, 100);
GAMMA = 1.4;

% Preallocating
v_f_2 = zeros(N,length(M_inf),length(AoA_2));
v_x_2 = zeros(N,length(M_inf),length(AoA_2));
v_z_2 = zeros(N,length(M_inf),length(AoA_2));
cp_2 = zeros(N,length(M_inf),length(AoA_2));
cl_2 = zeros(length(M_inf),length(AoA_2));
cm_0_2 = zeros(length(M_inf),length(AoA_2));
gamma_2 = zeros(N,length(M_inf),length(AoA_2));

for i = 1:length(AoA_2)
    for j = 1:length(M_inf)
        %Q_inf = 343*M_inf(j);
        Q_inf = 1;
        [v_f,v_x,v_z,cp,cl,cm_0,gamma] = Vortex(Q_inf,AoA_2(i),cosinus,sinus,longitud_panel,nodal_points,control_points,vector_tangent,c);
        v_f_2(:,j,i) = v_f;
        v_x_2(:,j,i) = v_x;
        v_z_2(:,j,i) = v_z;
        cp_2(:,j,i) = cp;
        cl_2(j,i) = cl;
        cm_0_2(j,i) = cm_0;
        gamma_2(:,j,i) = gamma;

        cp_min(i) = min(cp_2(:,j,i));

        cp_star(j,i) = (2/(GAMMA*M_inf(j).^2))*(((2+(GAMMA-1)*M_inf(j).^2)/(1+GAMMA))^(GAMMA/(GAMMA-1))-1);
        cp_Laitone(j,i) = cp_min(i)/(sqrt(1-M_inf(j)^2)+(cp_min(i)/2)*(M_inf(j)^2/sqrt(1-M_inf(j)^2))*(1+(GAMMA-1)*M_inf(j)^2/2));

        if abs(cp_Laitone(j,i)-cp_star(j,i))<0.2
            M_inf_2(i) = M_inf(j);
            break;
        end
    end
end
clear i;
%% APARTAT 3
AoA_3 = 6;
%M_inf_3 = [M_inf_final(2)-0.15 M_inf_final(2)-0.1 M_inf_final(2)-0.05 M_inf_final(2)];
M_inf_3 = linspace(0.01, 0.99, 100);
% Preallocating
v_f_3 = zeros(length(M_inf_3),N);
v_x_3 = zeros(length(M_inf_3),N);
v_z_3 = zeros(length(M_inf_3),N);
cp_3 = zeros(length(M_inf_3),N);
cl_3_incompressible = zeros(1,length(M_inf_3));
cm_0_3 = zeros(1,length(M_inf_3));
gamma_3 = zeros(length(M_inf_3),N);


for i = 1:length(M_inf_3)
    Q_inf=1;
    Beta(i) = sqrt(1-M_inf_3(i).^2);
    [v_f,v_x,v_z,cp,cl,cm_0,gamma] = Vortex(Q_inf,AoA_3,cosinus,sinus,longitud_panel,nodal_points,control_points,vector_tangent,c);
    v_f_3(i,:) = v_f;
    v_x_3(i,:) = v_x;
    v_z_3(i,:) = v_z;
    cp_3(i,:) = cp;
    cl_3_incompressible(i) = cl;
    cm_0_3(i) = cm_0;
    gamma_3(i,:) = gamma;
    cl_3(i) = cl_3_incompressible(i)./Beta(i);
end

%%
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
