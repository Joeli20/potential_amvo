function V = vortex_line(x,z,x_c,z_c,gamma,N)
R1_x = zeros(N-1,N-1);
R1_z = zeros(N-1,N-1);
R1_x = zeros(N-1,N-1);
R2_z = zeros(N-1,N-1);
    for i = 1:size(x_c,1)
        for j = 1:size(x_c,1)
            if i == j
                V(i,j,:) = 0;
            else
            R1_x(i,j) = x_c(i,1)-x(1,j);
            R1_z(i,j) = z_c(i,1)-z(1,j);
            R1(1,1) = R1_x(i,j);
            R1(1,3) = R1_z(i,j);
            R1(1,2) = 0;
            R2_x(i,j) = x_c(i,1)-x(1,j+1);
            R2_z(i,j) = z_c(i,1)-z(1,j+1);
            R2(1,1) = R2_x(i,j);
            R2(1,3) = R2_z(i,j);
            R2(1,2) = 0;
            R1_m(i,j) = sqrt((R1_x(i,j))^2+(R1_z(i,j))^2);
            R2_m(i,j) = sqrt((R2_x(i,j))^2+(R2_z(i,j))^2);
            L = cross(R1(1,:),R2(1,:))
            V(i,j,:) = (gamma(i)/(4*pi))*((R1_m(i,j)+R2_m(i,j))/R1_m(i,j))*R2_m(i,j)*(R1_m(i,j)*R2_m(i,j)+(R1_x(i,j)*R2_x(i,j)+R1_z(i,j)*R2_z(x,j)))*cross(R1,R2);
            end
        end
   end
end