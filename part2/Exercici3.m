clc; clear;
N = 64;
b=4.5;
rho = 1.225;
y = linspace(-2.25,2.25,N+1);
x = linspace(0,0,N+1);
z = linspace(1,1,N+1);
AoA_d = 0;
CL_0 = 0.24;
CL_alpha = 6.7;
x_c = zeros(N,1);
y_c = zeros(N,1);
z_c = ones(N,1);
gamma = zeros(N+1,1);
c = zeros(N,1);
theta = zeros(N,1);
c_t = 0.6;
c_r = 0.9;
theta_r = 0;
theta_t = 0;
m_chord_right = (c_t-c_r)/y(N+1);
n_chord_right = c_r;
m_theta_right = (theta_t-theta_r)/y(N+1);
n_theta_right = theta_r;
m_chord_left = (c_t-c_r)/y(1);
n_chord_left = c_r;
m_theta_left = (theta_t-theta_r)/y(1);
n_theta_left = theta_r;
AoA = AoA_d*(pi/180);
%Vel_inf = linspace(1,50,500);
Residu = zeros(500,1);
Res_minim = 50000;
for i = 1:N
    y_c(i) = (y(i)+y(i+1))/2;
    x_c(i) = (x(i)+x(i+1))/2;
    z_c(i) = (z(i)+z(i+1))/2;
    if i<(N+1)/2
    c(i) = m_chord_left*y_c(i)+n_chord_left;
    theta(i) = m_theta_left*y_c(i)+n_theta_left;
    else
    c(i) = m_chord_right*y_c(i)+n_chord_right;
    theta(i) = m_theta_right*y_c(i)+n_theta_right;
    end
end
S = c_t*b+((c_r-c_t)*(b))/2;
[V_inf1,V_inf2] = inf_vortex_line(x,y,z,x_c,y_c,z_c,N,AoA);
[V_vortex_line] = vortex_line(x,y,z,x_c,y_c,z_c,N);
%for j = 1:length(Vel_inf)
 %   Q_inf(1,1) = Vel_inf(1,j)*cos(AoA);
    Q_inf(1,1) = 33.02*cos(AoA);
    Q_inf(1,2) = 0;
    Q_inf(1,3) = 33.02*sin(AoA);
    %Q_inf(1,3) = Vel_inf(1,j)*sin(AoA);
    [gamma, V] = gamma_horsehoe(x_c,c,V_inf1,V_inf2,V_vortex_line,N,AoA,Q_inf,CL_0,CL_alpha,theta);
    for i = 1:N
        L(i,1) = gamma(i,1)*abs(y(1,i)-y(1,i+1)); %rho*norm(Q_inf)*
        C_l(i,1) = 2*((gamma(i,1)*abs(y(1,i)-y(1,i+1)))/(norm(Q_inf)*S));
        alpha_ind(i,1) = ((C_l(i,1)-CL_0)/CL_alpha)-AoA-theta(i,1);
        D_ind(i,1) = -rho*norm(Q_inf)*gamma(i,1)*abs(y(1,i)-y(1,i+1))*alpha_ind(i,1);
    end
    %Lift(j,1) = rho*norm(Q_inf)*sum(L);
    %Coef_Lift(j,1) = sum(C_l);
    %alpha_induit(j,1) = sum(alpha_ind);
    %Drag_ind(j,1) = sum(D_ind);
    Lift = rho*norm(Q_inf)*sum(L);
    Coef_Lift = sum(C_l);
    alpha_induit = sum(alpha_ind);
    Drag_ind = sum(D_ind);
    %Residu(j,1) = abs((40*9.81)-Lift(j,1));
    %if Residu(j,1) < Res_minim
     %   Res_minim = Residu(j,1);
      %  V_ok = Vel_inf(1,j);
    %end
%end
