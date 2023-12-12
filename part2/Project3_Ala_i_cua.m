clc;
clear;
%% Geometry
b=6;
c_root = 1.3;
c_tip = 0.7;
b_a = 2; %Span aileron, only one side
b_h = 2.2; %Span tail
c_root_h = 0.65;
c_tip_h = 0.45;
l = 3; %Distance wing vs tail
%NACA 0010 for wing and tail, incidence angle wing = 0;
i_ah = -2; % incidence angle tail in degrees;
%NACA 0015 for ailerons, and elevator.
Cm_025 = 0;
CL_0 = 0.24;
CL_alpha = 6.7;
CL_delta = 6.7;
rho = 1.225;
N = 128;
N_ala = N;
N_cua = N/2;
N_alerons = round((N_ala/2)*3/2);
Trigger = 0;
%% Discretization for normal use
x = linspace(0,0,N_ala+1);
y = linspace(-b/2,b/2,N_ala+1);
z = linspace(0,0,N_ala+1);
x_h = linspace(3,3,N_cua+1);
y_h = linspace(-b_h/2,b_h/2,N_cua+1);
z_h = linspace(0.05,0.05,N_cua+1);
x_c = zeros(N_ala,1);
y_c = zeros(N_ala,1);
z_c = zeros(N_ala,1);
x_c_h = zeros(N_cua,1);
y_c_h = zeros(N_cua,1);
z_c_h = zeros(N_cua,1);
c = zeros(N,1);
c_h = zeros(N_cua,1);
m_c_right = (c_tip-c_root)/y(N_ala+1);
n_c_right = c_root;
m_c_left = (c_tip-c_root)/y(1);
n_c_left = c_root;
m_c_h_right = (c_tip_h-c_root_h)/y_h(N_cua+1);
n_c_h_right = c_root_h;
m_c_h_left = (c_tip_h-c_root_h)/y_h(1);
n_c_h_left = c_root_h;
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
for i = 1:N_cua
    y_c_h(i) = (y_h(i)+y_h(i+1))/2;
    x_c_h(i) = (x_h(i)+x_h(i+1))/2;
    z_c_h(i) = (z_h(i)+z_h(i+1))/2;
    if i<(N_cua+1)/2
    c_h(i) = m_c_h_left*y_c_h(i)+n_c_h_left;
    else
    c_h(i) = m_c_h_right*y_c_h(i)+n_c_h_right;
    end
end
S_wing = (b*c_tip+(b*(c_root-c_tip)/2));
S_tail = (b_h*c_tip_h+(b_h*(c_root_h-c_tip_h)/2));
lambda_ala = c_tip/c_root;
lambda_tail = c_tip_h/c_root_h;
M_chord_ala = (2/3)*c_root*((1+lambda_ala+lambda_ala^2)/(1+lambda_ala));
M_chord_tail = (2/3)*c_root_h*((1+lambda_tail+lambda_tail^2)/(1+lambda_tail));

%% Wing and tail
xw_1(1,:) = c*3/4;
xw_2(1,:) = -c*1/4;
for i = 1:N_ala
    yw_b(i) = -b/2 + (i-1)*b/N_ala + (b/N_ala/2);
end
xt_1(1,:) = c_h*3/4;
xt_2(1,:) = -c_h*1/4;
for i = 1:N_cua
    yt_b(i) = -b_h/2 + (i-1)*b_h/N_cua +(b_h/N_cua/2);
end

%Wing
xw = [-c_tip*1/4 c_tip*3/4 xw_1 c_tip*3/4 -c_tip*1/4 ];
yw = [-b/2 -b/2 yw_b b/2 b/2 ];
figure
plot(xw,yw,Color='black')
hold on
plot(xw_2,yw_b,Color='black')
for i=1:N_ala
    plot([c(i)*3/4 -c(i)*1/4],[yw_b(i) yw_b(i)],Color='black')
end
%Tail
xt = [-c_tip_h*1/4 c_tip_h*3/4 xt_1 c_tip_h*3/4 -c_tip_h*1/4]+l;
yt = [-b_h/2 -b_h/2 yt_b b_h/2 b_h/2];
plot(xt,yt,Color='black')
hold on
plot(l+xt_2,yt_b,Color='black')
for i=1:N_cua
    plot(l+[c_h(i)*3/4 -c_h(i)*1/4],[yt_b(i) yt_b(i)],Color='black')
end
xlim([-1 4])

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

%[L_total,L_tail,Cl_c,L_wing,CL_wing, CL_tail, Cl_c_tail, CD_v_wing, cd_wing, cd_tail, CD_v_tail,CD_i_wing, CD_i_tail, CD_wing,CD_tail, Cd_i_wing,Cd_i_tail, Cd_v_wing,Cd_v_tail, alpha_i_wing,alpha_i_tail,Cm_0_ala,Cm_0_tail,Effi_ala,Effi_tail,Cm_roll_ala,Cm_roll_cua,M_roll_ala,M_roll_cua,N_pitch_cua,Cn_pitch_cua] = coefficients(rho,y,y_h,x_c,y_c,z_c,x_c_h,y_c_h,z_c_h,c,c_h,Q_inf,AoA, AoA_t,N_ala,N_cua,CL_0,CL_alpha,gamma,b,b_h, S_wing, S_tail,theta,Cd_v0,Cd_v_Cl,M_chord_ala,M_chord_tail, Cm_025);
%for i = 1:N_ala
%gamma_a(i,1) = gamma(i,1);
%end
%for i = N_ala+1:N_ala+N_cua
%gamma_t(i-N_ala,1) = gamma(i,1);
%end
%plot(y_c,gamma_a)
%hold on
%plot(y_c_h,gamma_t)

%% PART 1
theta = linspace(-6,6,10);
    for i=1:length(theta)
        theta_t(i,1) = deg2rad(theta(1,i));
        %Linear twist Calculation only in the wing
        theta_r = 0; %The wing in the middle has no twist.
        m_theta_right = (theta_t(i,1)-theta_r)/y(N);
        n_theta_right = theta_r;
        m_theta_left = (theta_t(i,1)-theta_r)/y(1);
        n_theta_left = theta_r;
        for j = 1:N_ala
            if j<(N+1)/2
                theta_dis(j,1) = m_theta_left*y_c(j,1)+n_theta_left;
            else
                theta_dis(j,1) = m_theta_right*y_c(j,1)+n_theta_right;
            end

        end
        Cd_v0 = 0.0075;
        Cd_v_Cl = 0.0055;
        delta_l = 0;
        delta_r = 0;
        delta_t = 0;
        AoA_d = 4;
        AoA_d_t = 4-i_ah;
        AoA = AoA_d*(pi/180);
        AoA_t = AoA_d_t*(pi/180);
        Q_inf(1,1) = 1*cos(AoA);
        Q_inf(1,2) = 0;
        Q_inf(1,3) = 1*sin(AoA);
        Q_inf_tail(1,1) = 1*cos(AoA_t);
        Q_inf_tail(1,2) = 0;
        Q_inf_tail(1,3) = 1*sin(AoA_t);

        %Wing-Wing
        [V_inf_1_1, V_inf_2_1] = inf_vortex_line(x,y,z,x_c,y_c,z_c,N_ala,N_ala,AoA);
        [V_AB_1] = vortex_line(x,y,z,x_c,y_c,z_c,N_ala,N_ala);
        %Wing-Tail
        [V_inf_1_2, V_inf_2_2] = inf_vortex_line(x_h,y_h,z_h,x_c,y_c,z_c,N_ala,N_cua,AoA);
        [V_AB_2] = vortex_line(x_h,y_h,z_h,x_c,y_c,z_c,N_ala,N_cua);
        %Tail-Wing
        [V_inf_1_3, V_inf_2_3] = inf_vortex_line(x,y,z,x_c_h,y_c_h,z_c_h,N_cua,N_ala,AoA_t);
        [V_AB_3] = vortex_line(x,y,z,x_c_h,y_c_h,z_c_h,N_cua,N_ala);
        %Tail-Tail
        [V_inf_1_4, V_inf_2_4] = inf_vortex_line(x_h,y_h,z_h,x_c_h,y_c_h,z_c_h,N_cua,N_cua,AoA_t);
        [V_AB_4] = vortex_line(x_h,y_h,z_h,x_c_h,y_c_h,z_c_h,N_cua,N_cua);


        [gamma(:,i), V_ij] = gamma_horsehoe_2 (c,c_h, V_inf_1_1,V_inf_2_1, V_inf_1_2,V_inf_2_2, V_inf_1_3,V_inf_2_3, V_inf_1_4,V_inf_2_4,V_AB_1,V_AB_2,V_AB_3,V_AB_4,N_ala,N_cua,N_alerons,AoA,AoA_t,Q_inf,CL_0,CL_alpha, CL_delta,theta_dis,delta_l,delta_r,delta_t);
        
        [L_total(i,1),L_tail(i,1),Cl_c,L_wing(i,1),CL_wing, CL_tail, Cl_c_tail, CD_v_wing, cd_wing, cd_tail, CD_v_tail,CD_i_wing, CD_i_tail, CD_wing,CD_tail, Cd_i_wing,Cd_i_tail, Cd_v_wing,Cd_v_tail, alpha_i_wing,alpha_i_tail,Cm_0_ala,Cm_0_tail,Effi_ala(i,1),Effi_tail(i,1),Cm_roll_ala,Cm_roll_cua,M_roll_ala,M_roll_cua,N_pitch_cua,Cn_pitch_cua] = coefficients(rho,y,y_h,x_c,y_c,z_c,x_c_h,y_c_h,z_c_h,c,c_h,Q_inf,AoA, AoA_t,N_ala,N_cua,CL_0,CL_alpha,gamma,b,b_h, S_wing, S_tail,theta_dis,Cd_v0,Cd_v_Cl,M_chord_ala,M_chord_tail, Cm_025);
        if Effi_ala(i,1)>Trigger
            Trigger = Effi_ala(i,1);
            THETA = theta_dis;
            THETA_T = theta_t(i,1);
            THETA_T_d = rad2deg(THETA_T);
            Cl_c_bo = Cl_c;
            L_total_bo = L_total(i,1);
        end
    end