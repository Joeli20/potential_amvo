function V_ab = vortex_line(x,y,z,x_c,y_c,z_c,N_1,N_2)
V_ab = zeros(N_1,N_2,3);
    for i = 1:N_1
        for j = 1:N_2
            R1 = [x_c(i)-x(j)   y_c(i)-y(j) z_c(i)-z(j)];
            R2 = [x_c(i)-x(j+1)   y_c(i)-y(j+1) z_c(i)-z(j+1)];
            R1_n = norm(R1);
            R2_n = norm(R2);
            if (norm(cross(R1,R2)) < 0.00001)
                V_ab(i,j,:) = [0 0 0];
            else
                V_ab(i,j,:) = 1/(4*pi)*(R1_n+R2_n)/(R1_n*R2_n*(R1_n*R2_n+dot(R1,R2)))*cross(R1,R2); 
            end
        end
   end
end