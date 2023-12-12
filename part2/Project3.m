clc;
clear;
close all;

%% INPUT DATA
% Geometry
b = 6; % Wing Span
c_r = 1.3; % Root Chord
c_t = 0.7; % Tip Chord
b_a = 2; % Aileron Span, only one side
b_h = 2.2; % Tail Span
c_r_h = 0.65; % Tail Root Chord
c_t_h = 0.45; % Tail Tip Chord
l = 3; % Distance Wing - tail

%Airfoil parameters
% NACA 0010 for wing and tail, incidence angle wing = 0;
% NACA 0015 for ailerons, and elevator.

i_w = 0; % Incidence angle Wing in degrees;
i_h = -2; % Incidence angle tail in degrees;

Cm_025 = 0; % NACA 0010
Cl_0 = 0; % NACA 0010
Cl_alpha = 6.7265; % NACA 0010
Cl_delta = 4.1554; % NACA 0015

N = 256; % Sections INPUT
N_w = N; % Wing sections
N_h = N/2; % Tail sections
N_a = round((N_w/2)*3/2); % Number of ailerons

% Other fluid parameters
rho = 1.225;

%% DISCRETITZATION
% Preallocating
x = linspace(0,0,N_w+1);
y = linspace(-b/2,b/2,N_w+1);
z = linspace(0,0,N_w+1);
x_h = linspace(3,3,N_h+1);
y_h = linspace(-b_h/2,b_h/2,N_h+1);
z_h = linspace(-0.005,-0.005,N_h+1);
x_c = zeros(N_w,1);
y_c = zeros(N_w,1);
z_c = zeros(N_w,1);
x_c_h = zeros(N_h,1);
y_c_h = zeros(N_h,1);
z_c_h = zeros(N_h,1);
c = zeros(N,1);
c_h = zeros(N_h,1);

% Previous calculations
m_c_right = (c_t-c_r)/y(N_w+1);
n_c_right = c_r;
m_c_left = (c_t-c_r)/y(1);
n_c_left = c_r;
m_c_h_right = (c_t_h-c_r_h)/y_h(N_h+1);
n_c_h_right = c_r_h;
m_c_h_left = (c_t_h-c_r_h)/y_h(1);
n_c_h_left = c_r_h;

for i = 1:N
    y_c(i) = (y(i)+y(i+1))/2;
    x_c(i) = (x(i)+x(i+1))/2;
    z_c(i) = (z(i)+z(i+1))/2;
    if i<(N+1)/2
    c(i) = m_c_left*y_c(i)+n_c_left;
    else
    c(i) = m_c_right*y_c(i)+n_c_right;
    end
end
for i = 1:N_h
    y_c_h(i) = (y_h(i)+y_h(i+1))/2;
    x_c_h(i) = (x_h(i)+x_h(i+1))/2;
    z_c_h(i) = (z_h(i)+z_h(i+1))/2;
    if i<(N_h+1)/2
    c_h(i) = m_c_h_left*y_c_h(i)+n_c_h_left;
    else
    c_h(i) = m_c_h_right*y_c_h(i)+n_c_h_right;
    end
end

% Post calculations
S_w = (b*c_t+(b*(c_r-c_t)/2));
S_h = (b_h*c_t_h+(b_h*(c_r_h-c_t_h)/2));
lambda_w = c_t/c_r;
lambda_h = c_t_h/c_r_h;
m_chord_w = (2/3)*c_r*((1+lambda_w+lambda_w^2)/(1+lambda_w));
m_chord_h = (2/3)*c_r_h*((1+lambda_h+lambda_h^2)/(1+lambda_h));

%% PART 1

% Preallocating
theta = linspace(-5,5,10);
twist = theta;
theta_t = zeros(length(theta),1);
theta_dis = zeros(N_w,1);
L_total = zeros(length(theta),1);
L_w = zeros(length(theta),1);
L_h = zeros(length(theta),1);
e_w = zeros(length(theta),1);
e_h = zeros(length(theta),1);

Trigger = 0;

for i=1:length(theta)
    theta_t(i,1) = deg2rad(theta(1,i));
    %Linear twist Calculation only in the wing
    theta_r = 0; % No twist on wing center
    m_theta_right = (theta_t(i,1)-theta_r)/y(N);
    n_theta_right = theta_r;
    m_theta_left = (theta_t(i,1)-theta_r)/y(1);
    n_theta_left = theta_r;
    for j = 1:N_w
        if j<(N+1)/2
            theta_dis(j,1) = m_theta_left*y_c(j,1)+n_theta_left;
        else
            theta_dis(j,1) = m_theta_right*y_c(j,1)+n_theta_right;
        end

    end
    
    % Aerodynamic parameters
    Cd_0 = 0.0075;
    Cd_Cl = 0.0055;
    delta_l = 0;
    delta_r = 0;
    delta_t = 0;
    AoA_d = 4;
    AoA_d_t = 4+i_h;
    AoA = AoA_d*(pi/180);
    AoA_t = AoA_d_t*(pi/180);
    Q_inf(1,1) = 1*cos(AoA);
    Q_inf(1,2) = 0;
    Q_inf(1,3) = 1*sin(AoA);
    Q_inf_h(1,1) = 1*cos(AoA_t);
    Q_inf_h(1,2) = 0;
    Q_inf_h(1,3) = 1*sin(AoA_t);

    % INTERACTIONS
    % Wing-Wing
    [V_inf_1_1, V_inf_2_1] = inf_vortex_line(x,y,z,x_c,y_c,z_c,N_w,N_w,AoA);
    [V_AB_1] = vortex_line(x,y,z,x_c,y_c,z_c,N_w,N_w);
    % Wing-Tail
    [V_inf_1_2, V_inf_2_2] = inf_vortex_line(x_h,y_h,z_h,x_c,y_c,z_c,N_w,N_h,AoA);
    [V_AB_2] = vortex_line(x_h,y_h,z_h,x_c,y_c,z_c,N_w,N_h);
    % Tail-Wing
    [V_inf_1_3, V_inf_2_3] = inf_vortex_line(x,y,z,x_c_h,y_c_h,z_c_h,N_h,N_w,AoA);
    [V_AB_3] = vortex_line(x,y,z,x_c_h,y_c_h,z_c_h,N_h,N_w);
    % Tail-Tail
    [V_inf_1_4, V_inf_2_4] = inf_vortex_line(x_h,y_h,z_h,x_c_h,y_c_h,z_c_h,N_h,N_h,AoA);
    [V_AB_4] = vortex_line(x_h,y_h,z_h,x_c_h,y_c_h,z_c_h,N_h,N_h);

    % Gamma calculation
    [gamma(:,i), V_ij] = gamma_horsehoe_2 (c,c_h, V_inf_1_1,V_inf_2_1, V_inf_1_2,V_inf_2_2, V_inf_1_3, ...
        V_inf_2_3, V_inf_1_4,V_inf_2_4,V_AB_1,V_AB_2,V_AB_3,V_AB_4,N_w,N_h,N_a,AoA,AoA_t,Q_inf,Cl_0, ...
        Cl_alpha, Cl_delta,theta_dis,delta_l,delta_r,delta_t);
    
    % Coefficients calculation
    [coef] = coefficients(rho,y,y_h,x_c,y_c, ...
        z_c,x_c_h,y_c_h,z_c_h,c,c_h,Q_inf,AoA, AoA_t,N_w,N_h,Cl_0,Cl_alpha,gamma(:,i),b,b_h, S_w, S_h, ...
        theta_dis,Cd_0,Cd_Cl,m_chord_w,m_chord_h, Cm_025);

    L_total(i,1) = coef.L_total;
    e_w(i,1) = coef.Effi_ala;
    Cl_c(:,i) = coef.Cl_c_wing;
    Cl_c_tail(:,i) = coef.Cl_c_tail;
    Cd_c_v(:,i) = coef.Cd_v_wing;
    Cd_c_i(:,i) = coef.Cd_i_wing;
    Cd_c_v_tail(:,i) = coef.Cd_v_tail;
    Cd_c_i_tail(:,i) = coef.Cd_i_tail;
    Alpha_c_i(:,i) = coef.alpha_i_wing;
    Alpha_c_i_tail(:,i) = coef.alpha_i_tail;
    CL(:,i) = coef.CL_wing;
    CL_tail(:,i) = coef.CL_tail;
    CD(:,i) = coef.CD_wing;
    CD_tail(:,i) = coef.CD_tail;
    CM(:,i) = coef.Cm_0_ala;
    CM_tail(:,i) = coef.Cm_0_tail;

    if e_w(i,1)>Trigger
        Trigger = e_w(i,1);
        THETA = theta_dis;
        THETA_T = theta_t(i,1);
        THETA_T_d = rad2deg(THETA_T);
        Cl_c_bo = Cl_c;
        L = L_total(i,1);
        k = i;
    end
end

plot(y_c/(b/2),Cl_c(:,k),y_c_h/(b/2),Cl_c_tail(:,k))
grid on
plot(y_c/(b/2),Cd_c_v(:,k),y_c_h/(b/2),Cd_c_v_tail(:,k));
plot(y_c/(b/2),Cd_c_i(:,k),y_c_h/(b/2),Cd_c_i_tail(:,k));
plot(y_c/(b/2),Alpha_c_i(:,k),y_c_h/(b/2),Alpha_c_i_tail(:,k));

%% Part 2

    % Aerodynamic parameters
    Cd_0 = 0.0075;
    Cd_Cl = 0.0055;
    delta_l = 0;
    delta_r = 0;
    delta_t = deg2rad(15);
    AoA_d = 4;
    AoA_d_t = 4+i_h;
    AoA = AoA_d*(pi/180);
    AoA_t = AoA_d_t*(pi/180);
    Q_inf(1,1) = 1*cos(AoA);
    Q_inf(1,2) = 0;
    Q_inf(1,3) = 1*sin(AoA);
    Q_inf_h(1,1) = 1*cos(AoA_t);
    Q_inf_h(1,2) = 0;
    Q_inf_h(1,3) = 1*sin(AoA_t);
    theta = zeros(N_w,1);


    % INTERACTIONS
    % Wing-Wing
    [V_inf_1_1, V_inf_2_1] = inf_vortex_line(x,y,z,x_c,y_c,z_c,N_w,N_w,AoA);
    [V_AB_1] = vortex_line(x,y,z,x_c,y_c,z_c,N_w,N_w);
    % Wing-Tail
    [V_inf_1_2, V_inf_2_2] = inf_vortex_line(x_h,y_h,z_h,x_c,y_c,z_c,N_w,N_h,AoA);
    [V_AB_2] = vortex_line(x_h,y_h,z_h,x_c,y_c,z_c,N_w,N_h);
    % Tail-Wing
    [V_inf_1_3, V_inf_2_3] = inf_vortex_line(x,y,z,x_c_h,y_c_h,z_c_h,N_h,N_w,AoA);
    [V_AB_3] = vortex_line(x,y,z,x_c_h,y_c_h,z_c_h,N_h,N_w);
    % Tail-Tail
    [V_inf_1_4, V_inf_2_4] = inf_vortex_line(x_h,y_h,z_h,x_c_h,y_c_h,z_c_h,N_h,N_h,AoA);
    [V_AB_4] = vortex_line(x_h,y_h,z_h,x_c_h,y_c_h,z_c_h,N_h,N_h);

    % Gamma calculation
    [gamma_2, V_ij_2] = gamma_horsehoe_2 (c,c_h, V_inf_1_1,V_inf_2_1, V_inf_1_2,V_inf_2_2, V_inf_1_3, ...
        V_inf_2_3, V_inf_1_4,V_inf_2_4,V_AB_1,V_AB_2,V_AB_3,V_AB_4,N_w,N_h,N_a,AoA,AoA_t,Q_inf,Cl_0, ...
        Cl_alpha, Cl_delta,theta_dis,delta_l,delta_r,delta_t);
    
    % Coefficients calculation
    [coef_2] = coefficients(rho,y,y_h,x_c,y_c, ...
        z_c,x_c_h,y_c_h,z_c_h,c,c_h,Q_inf,AoA, AoA_t,N_w,N_h,Cl_0,Cl_alpha,gamma_2,b,b_h, S_w, S_h, ...
        theta_dis,Cd_0,Cd_Cl,m_chord_w,m_chord_h, Cm_025);

    Cl_c_2(:,i) = coef_2.Cl_c_wing;
    Cl_c_tail_2(:,i) = coef_2.Cl_c_tail;

%% Part 3

    % Aerodynamic parameters
    Cd_0 = 0.0075;
    Cd_Cl = 0.0055;
    delta_l = deg2rad(-10);
    delta_r = deg2rad(10);
    delta_t = deg2rad(0);
    AoA_d = 4;
    AoA_d_t = 4+i_h;
    AoA = AoA_d*(pi/180);
    AoA_t = AoA_d_t*(pi/180);
    Q_inf(1,1) = 1*cos(AoA);
    Q_inf(1,2) = 0;
    Q_inf(1,3) = 1*sin(AoA);
    Q_inf_h(1,1) = 1*cos(AoA_t);
    Q_inf_h(1,2) = 0;
    Q_inf_h(1,3) = 1*sin(AoA_t);
    theta = zeros(N_w,1);


    % INTERACTIONS
    % Wing-Wing
    [V_inf_1_1, V_inf_2_1] = inf_vortex_line(x,y,z,x_c,y_c,z_c,N_w,N_w,AoA);
    [V_AB_1] = vortex_line(x,y,z,x_c,y_c,z_c,N_w,N_w);
    % Wing-Tail
    [V_inf_1_2, V_inf_2_2] = inf_vortex_line(x_h,y_h,z_h,x_c,y_c,z_c,N_w,N_h,AoA);
    [V_AB_2] = vortex_line(x_h,y_h,z_h,x_c,y_c,z_c,N_w,N_h);
    % Tail-Wing
    [V_inf_1_3, V_inf_2_3] = inf_vortex_line(x,y,z,x_c_h,y_c_h,z_c_h,N_h,N_w,AoA);
    [V_AB_3] = vortex_line(x,y,z,x_c_h,y_c_h,z_c_h,N_h,N_w);
    % Tail-Tail
    [V_inf_1_4, V_inf_2_4] = inf_vortex_line(x_h,y_h,z_h,x_c_h,y_c_h,z_c_h,N_h,N_h,AoA);
    [V_AB_4] = vortex_line(x_h,y_h,z_h,x_c_h,y_c_h,z_c_h,N_h,N_h);

    % Gamma calculation
    [gamma_3, V_ij_3] = gamma_horsehoe_2 (c,c_h, V_inf_1_1,V_inf_2_1, V_inf_1_2,V_inf_2_2, V_inf_1_3, ...
        V_inf_2_3, V_inf_1_4,V_inf_2_4,V_AB_1,V_AB_2,V_AB_3,V_AB_4,N_w,N_h,N_a,AoA,AoA_t,Q_inf,Cl_0, ...
        Cl_alpha, Cl_delta,theta_dis,delta_l,delta_r,delta_t);
    
    % Coefficients calculation
    [coef_3] = coefficients(rho,y,y_h,x_c,y_c, ...
        z_c,x_c_h,y_c_h,z_c_h,c,c_h,Q_inf,AoA, AoA_t,N_w,N_h,Cl_0,Cl_alpha,gamma_3,b,b_h, S_w, S_h, ...
        theta_dis,Cd_0,Cd_Cl,m_chord_w,m_chord_h, Cm_025);

    Cl_c_3(:,i) = coef_3.Cl_c_wing;
    Cl_c_tail_3(:,i) = coef_3.Cl_c_tail;

%% CODE END
%% TEST

%Cd_v0 = 0.0075;
%Cd_v_Cl = 0.0055;
%AoA_d = 4;
%AoA_d_t = 4-i_ah;
%AoA = AoA_d*(pi/180);
%AoA_t = AoA_d_t*(pi/180);
%Q_inf(1,1) = 1*cos(AoA);
%Q_inf(1,2) = 0;
%Q_inf(1,3) = 1*sin(AoA);
%Q_inf_tail(1,1) = 1*cos(AoA_t);
%Q_inf_tail(1,2) = 0;
%Q_inf_tail(1,3) = 1*sin(AoA_t);
%theta = zeros(N_ala+N_cua,1);

%Wing-Wing
%[V_inf_1_1, V_inf_2_1] = inf_vortex_line(x,y,z,x_c,y_c,z_c,N_ala,N_ala,AoA);
%[V_AB_1] = vortex_line(x,y,z,x_c,y_c,z_c,N_ala,N_ala);
%Wing-Tail
%[V_inf_1_2, V_inf_2_2] = inf_vortex_line(x_h,y_h,z_h,x_c,y_c,z_c,N_ala,N_cua,AoA);
%[V_AB_2] = vortex_line(x_h,y_h,z_h,x_c,y_c,z_c,N_ala,N_cua);
%Tail-Wing
%[V_inf_1_3, V_inf_2_3] = inf_vortex_line(x,y,z,x_c_h,y_c_h,z_c_h,N_cua,N_ala,AoA_t);
%[V_AB_3] = vortex_line(x,y,z,x_c_h,y_c_h,z_c_h,N_cua,N_ala);
%Tail-Tail
%[V_inf_1_4, V_inf_2_4] = inf_vortex_line(x_h,y_h,z_h,x_c_h,y_c_h,z_c_h,N_cua,N_cua,AoA_t);
%[V_AB_4] = vortex_line(x_h,y_h,z_h,x_c_h,y_c_h,z_c_h,N_cua,N_cua);

%[gamma, V_ij] = gamma_horsehoe_2 (c,c_h, V_inf_1_1,V_inf_2_1, V_inf_1_2,V_inf_2_2, V_inf_1_3,V_inf_2_3, V_inf_1_4,V_inf_2_4,V_AB_1,V_AB_2,V_AB_3,V_AB_4,N_ala,N_cua,N_alerons,AoA,AoA_t,Q_inf,CL_0,CL_alpha, CL_delta);

%[coef] = coefficients(rho,y,y_h,x_c,y_c,z_c,x_c_h,y_c_h,z_c_h,c,c_h,Q_inf,AoA, AoA_t,N_ala,N_cua,CL_0,CL_alpha,gamma,b,b_h, S_wing, S_tail,theta,Cd_v0,Cd_v_Cl,M_chord_ala,M_chord_tail, Cm_025);
%for i = 1:N_ala
%gamma_a(i,1) = gamma(i,1);
%end
%for i = N_ala+1:N_ala+N_cua
%gamma_t(i-N_ala,1) = gamma(i,1);
%end
%plot(y_c,gamma_a)
%hold on
%plot(y_c_h,gamma_t)