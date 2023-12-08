N = 9;
y = linspace(-2.25,2.25,N);
x = linspace(0,0,N);
z = linspace(1,1,N);
Vel_inf = 10;
AoA_d = 0;
CL_0 = 0.24;
CL_alpha = 6.7;
x_c = zeros(N-1,1);
y_c = zeros(N-1,1);
z_c = ones(N-1,1);
gamma = ones(N-1,1);
c = zeros(N-1,1);
theta = zeros(N-1,1);
c_t = 0.6;
c_r = 0.9;
theta_r = 0;
theta_t = 0.3;
m_chord_right = (c_t-c_r)/y(N);
n_chord_right = c_r;
m_theta_right = (theta_t-theta_r)/y(N);
n_theta_right = theta_r;
m_chord_left = (c_t-c_r)/y(1);
n_chord_left = c_r;
m_theta_left = (theta_t-theta_r)/y(1);
n_theta_left = theta_r;
AoA = AoA_d*(pi/180);
Q_inf(1,1) = cos(AoA);
Q_inf(1,2) = 0;
Q_inf(1,3) = sin(AoA);
for i = 1:N-1
    y_c(i) = (y(i)+y(i+1))/2;
    x_c(i) = (x(i)+x(i+1))/2;
    z_c(i) = (z(i)+z(i+1))/2;
    if i<N/2
    c(i) = m_chord_left*y_c(i)+n_chord_left;
    theta(i) = m_theta_left*y_c(i)+n_theta_left;
    else
    c(i) = m_chord_right*y_c(i)+n_chord_right;
    theta(i) = m_theta_right*y_c(i)+n_theta_right;
    end
end

[V_inf1,V_inf2] = inf_vortex_line(x,y,z,x_c,y_c,z_c,N,gamma,AoA);
[V_vortex_line] = vortex_line(x,y,z,x_c,y_c,z_c,gamma,N);
[gamma, V] = gamma_horsehoe(x_c,c,V_inf1,V_inf2,V_vortex_line,N,AoA,Q_inf,CL_0,CL_alpha,theta);