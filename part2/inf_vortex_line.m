function [V_inf1,V_inf2] = inf_vortex_line(x,y,z,x_c,y_c,z_c,N,gamma,AoA)
%R1 = zeros(N-1,N-1,3);
%R2 = zeros(N-1,N-1,3);
%ur = zeros(N-1,N-1,3);
%ur_2 = zeros(N-1,N-1,3);
%ur_1 = zeros(N-1,N-1,3);
V_inf1 = zeros(N-1,N-1,3);
V_inf2 = zeros(N-1,N-1,3);
    for i = 1:size(y_c,1)
        for j = 1:size(y_c,1)
            %Seminfinite 1
            R2 = [x_c(i)-x(j)    y_c(i)-y(j) z_c(i)-z(j)];
            ur = [-cosd(AoA) 0   -sind(AoA)];
            ur_2 = R2./norm(R2);
            if (norm(cross(ur,R2)) < 0.00001)
            V_inf1(i,j,:)= [0 0 0];
            else
            V_inf1(i,j,:) = 1/(4*pi)*(1-dot(ur,ur_2))/(norm(cross(ur,R2))^2)*cross(ur,R2); 
            end
            %Seminfinite 2
            R1 = [x_c(i)-x(j+1)  y_c(i)-y(j+1)   z_c(i)-z(j+1)];
            ur_1 = R1./norm(R2);
            if (norm(cross(ur,R1)) < 0.00001)
            V_inf2(i,j,:)= [0 0 0];
            else
            V_inf2(i,j,:) = 1/(4*pi)*(1-dot(ur,ur_1))/(norm(cross(ur,R1))^2)*cross(ur,R1); 
            end
        end
   end
end