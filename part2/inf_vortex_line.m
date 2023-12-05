function V = inf_vortex_line(x,z,x_c,z_c,gamma)
    for i = 1:length(x_c,1)
        for j = 1:length(x_c,1)
            R1_x(i,j) = x_c(i,1)-x(j,1);
            R1_z(i,j) = z_c(i,1)-z(j,1);
            R1(1,1) = R1_x(i,j);
            R1(1,2) = R1_z(i,j);
            R2_x(i,j) = -x_c(i,1)+x(j+1,1);
            R2_z(i,j) = -z_c(i,1)+z(j+1,1);
            R2(1,1) = R2_x(i,j);
            R2(1,2) = R2_z(i,j);
            R1_m(i,j) = sqrt((R1_x(i,j))^2+(R1_z(i,j))^2);
            R2_m(i,j) = sqrt((R2_x(i,j))^2+(R2_z(i,j))^2);
            Ur1_x(i,j) = R1_x(i,j)/R1_m(i,j);
            Ur1_z(i,j) = R1_z(i,j)/R1_m(i,j);
            Ur1(1,1) = Ur1_x(i,j);
            Ur1(1,2) = Ur1_z(i,j);
            Ur2_x(i,j) = R2_x(i,j)/R1_m(i,j);
            Ur2_z(i,j) = R2_z(i,j)/R1_m(i,j);
            Ur2(1,1) = Ur2_x(i,j);
            Ur2(1,2) = Ur2_z(i,j);
            V(i,j,:) = (gamma(i)/(4*pi))*((1-(Ur1_x(i,j)*Ur2_x(i,j)+Ur1_z(i,j)*Ur2_z(i,j)))/(norm(cross(Ur2,R2)))^2)*cross(Ur1,R2);
        end
   end
end