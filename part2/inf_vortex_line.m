function [V_inf1,V_inf2] = inf_vortex_line(x,y,z,x_c,y_c,z_c,N_1,N_2,AoA)
% Calcula la influencia entre nodes d'una recta de vortex infinita
%
% Escrit per: Joel Campo, Jordi Gallart, Martí Santamaria, 2023
% Group 16. AMVO. MUEA.
%
V_inf1 = zeros(N_1,N_2,3);
V_inf2 = zeros(N_1,N_2,3);
    for i = 1:N_1
        for j = 1:N_2
            %Seminfinite 1
            R2 = [x_c(i,1)-x(1,j)    y_c(i,1)-y(1,j) z_c(i,1)-z(1,j)];
            ur = [-cos(AoA) 0   -sin(AoA)];
            ur_2 = R2./norm(R2);
            if (norm(cross(ur,R2)) < 0.00001)
            V_inf1(i,j,:)= [0 0 0];
            else
            V_inf1(i,j,:) = 1/(4*pi)*(1-dot(ur,ur_2))/(norm(cross(ur,R2))^2)*cross(ur,R2); 
            end
            %Seminfinite 2
            R1 = [x_c(i,1)-x(1,j+1)  y_c(i,1)-y(1,j+1)   z_c(i,1)-z(1,j+1)];
            ur_1 = R1./norm(R2);
            if (norm(cross(ur,R1)) < 0.00001)
            V_inf2(i,j,:)= [0 0 0];
            else
            V_inf2(i,j,:) = 1/(4*pi)*(1-dot(ur,ur_1))/(norm(cross(ur,R1))^2)*cross(ur,R1); 
            end
        end
   end
end