function [V_f,P_f] = ViP_fons(Q_inf,alpha,l_p,nodes,n)

a = zeros(size(n,2),size(n,2));
b = zeros(size(n,2),1);
sigma = zeros(size(n,2),1);

for i = 1:size(n,2)
    for j = 1:size(n,2)
        b(i,1) = -Q_inf*n(i);
        if i ~= j
            x_pan(1,j) = (x_ci-nodes(1,j))*cos(alpha(i))-(z_ci-nodes(2,j))*sin(alpha(i));
            z_pan(1,j) = (x_ci-nodes(1,j))*sin(alpha(i))-(z_ci-nodes(2,j))*cos(alpha(i));
            r1 = sqrt(x_pan^2+z_pan^2);
            r2 = sqrt((x_pan-l_p)^2+z_pan^2);
            u_ind_p = 1/(4*pi)*log((r1^2)/(r2^2)); 
            theta1 = atan(z_pan/x_pan);
            theta2 = atan(z_pan/(x_pan-l_p));
            w_ind_p = (theta2-theta1)/(2*pi);
            V = [u_ind_p; w_ind_p];
            a(i,j) = V(i,j)*n(i);
        end
        if i == j
            a(i,j) = 1/2;
        end
        sigma(i,1) = a(i,j)\b(i,1);
    end
end
V_f = Q_inf+sum(sigma(:,1)*V(:,1));
P_f = 1-(norm(V_f)/Q_inf)^2;

end