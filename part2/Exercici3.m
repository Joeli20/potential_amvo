N = 9;
x = linspace(-2.25,2.25,N);
x_c = zeros(N-1,1);
c = zeros(N-1,1);
theta = zeros(N-1,1);
c_t = 0.6;
c_r = 0.9;
theta_r = 0;
theta_t = 0.3;
m_chord_right = (c_t-c_r)/x(N);
n_chord_right = c_r;
m_theta_right = (theta_t-theta_r)/x(N);
n_theta_right = theta_r;
m_chord_left = (c_t-c_r)/x(1);
n_chord_left = c_r;
m_theta_left = (theta_t-theta_r)/x(1);
n_theta_left = theta_r;
for i = 1:N-1
    x_c(i) = (x(i)+x(i+1))/2;
    if i<N/2
    c(i) = m_chord_left*x_c(i)+n_chord_left;
    theta(i) = m_theta_left*x_c(i)+n_chord_left;
    else
    c(i) = m_chord_right*x_c(i)+n_chord_right;
    theta(i) = m_theta_right*x_c(i)+n_chord_right;
    end
end