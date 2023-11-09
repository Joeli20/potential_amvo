function [V_f,P_f,a, sigma, V] = ViP_fons(Q_inf,cosinus,sinus,l_p,nodes,center,normal)

a = zeros(size(normal,2),size(normal,2));
b = zeros(size(normal,2),1);
sigma = zeros(size(normal,2),1);
V = zeros(size(normal,2),1);

for i = 1:size(normal,2)
    for j = 1:size(normal,2)
        b(i,1) = -Q_inf*normal(i);
        if i ~= j
            x_pan(1,j) = (center(1,i)-nodes(1,j))*cosinus(1,j)-(center(2,i)-nodes(2,j))*sinus(1,j);
            z_pan(1,j) = (center(1,i)-nodes(1,j))*sinus(1,j)-(center(2,i)-nodes(2,j))*cosinus(1,j);
            r1(i,j) = sqrt(x_pan(1,j).^2+z_pan(1,j).^2);
            r2(i,j) = sqrt((x_pan(1,j)-l_p(j,1)).^2+z_pan(1,j).^2);
            u_ind_p(i,j) = 1/(4*pi)*log((r1(i,j).^2)/(r2(i,j).^2)); 
            theta1(i,j) = atan2(z_pan(1,j),x_pan(1,j));
            theta2(i,j) = atan2(z_pan(1,j),(x_pan(1,j)-l_p(j,1)));
            w_ind_p(i,j) = (theta2(i,j)-theta1(i,j))/(2*pi);
            u_ind_global(i,j) = u_ind_p(i,j)*cosinus(1,j) + w_ind_p(i,j)*sinus(1,j);
            w_ind_global(i,j) = -u_ind_p(i,j)*sinus(1,j) + w_ind_p(i,j)*cosinus(1,j);
            V(i,1) = sqrt(u_ind_global(i,j).^2+w_ind_global(i,j).^2);
            a(i,j) = u_ind_global(i,j)*normal(1,i)+w_ind_global(i,j)*normal(2,i);
        end
        if i == j
            a(i,j) = 1/2;
        end
        sigma(i,1) = a(i,j)\b(i,1);
    end
    V_f(i) = Q_inf+sum(sigma(:,1).*V(:,1));
    P_f(i) = 1-((V_f(i)).^2)./(Q_inf^2);
end


end