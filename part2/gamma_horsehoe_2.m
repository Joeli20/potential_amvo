function [gamma, V_ij] = gamma_horsehoe_2 (c,c_h, Vinf1_1,Vinf2_1, Vinf1_2,Vinf2_2, Vinf1_3,Vinf2_3, Vinf1_4,Vinf2_4,V_AB_1,V_AB_2,V_AB_3,V_AB_4,N_1,N_2,N_3,AoA,AoA_t,Q_inf,Cl_0,Cl_alpha,Cl_delta, theta, delta_l,delta_r,delta_t)
% Aplica la teoria de Horsehoe i després es calcula Gamma
%
% Escrit per: Joel Campo, Jordi Gallart, Martí Santamaria, 2023
% Group 16. AMVO. MUEA.
%
V_ij = zeros(N_1,N_2,3);
a = zeros(N_1+N_2,N_1+N_2);
b = zeros(N_1+N_2,1);
for i = 1:N_1
    for j = 1:N_1
            V_ij(i,j,:) = V_AB_1(i,j,:)+Vinf1_1(i,j,:)-Vinf2_1(i,j,:);
        if i == j
            a(i,j) = -0.5*Cl_alpha*c(i,1)*V_ij(i,j,3)+1;
        else
            a(i,j) = -0.5*Cl_alpha*c(i,1)*V_ij(i,j,3);
        end
    end
    for j = (N_1+1):(N_1+N_2)
        V_ij(i,j,:) = V_AB_2(i,j-N_1,:)+Vinf1_2(i,j-N_1,:)-Vinf2_2(i,j-N_1,:);
        if i == j
            a(i,j) = -0.5*Cl_alpha*c(i,1)*V_ij(i,j,3)+1;
        else
            a(i,j) = -0.5*Cl_alpha*c(i,1)*V_ij(i,j,3);
        end
    end
end
for i = (N_1+1):(N_1+N_2)
    for j = 1:N_1
        V_ij(i,j,:) = V_AB_3(i-N_1,j,:)+Vinf1_3(i-N_1,j,:)-Vinf2_3(i-N_1,j,:);
        if i == j
            a(i,j) = -0.5*Cl_alpha*c_h(i-N_1,1)*V_ij(i,j,3)+1;
        else
            a(i,j) = -0.5*Cl_alpha*c_h(i-N_1,1)*V_ij(i,j,3);
        end
    end
    for j = (N_1+1):(N_1+N_2)
        V_ij(i,j,:) = V_AB_4(i-N_1,j-N_1,:)+Vinf1_4(i-N_1,j-N_1,:)-Vinf2_4(i-N_1,j-N_1,:);
        if i == j
            a(i,j) = -0.5*Cl_alpha*c_h(i-N_1,1)*V_ij(i,j,3)+1;
        else
            a(i,j) = -0.5*Cl_alpha*c_h(i-N_1,1)*V_ij(i,j,3);
        end
    end
for i = 1:N_3
    b(i,1) = 0.5*c(i,1)*norm(Q_inf)*(Cl_0+Cl_alpha*(AoA+theta(i,1))+Cl_delta*delta_l);
end
for i = N_3+1:N_1-N_3
    b(i,1) = 0.5*c(i,1)*norm(Q_inf)*(Cl_0+Cl_alpha*(AoA+theta(i,1)));
end
for i = N_1-N_3+1:N_1
    b(i,1) = 0.5*c(i,1)*norm(Q_inf)*(Cl_0+Cl_alpha*(AoA+theta(i,1))+Cl_delta*delta_r);
end
for i = N_1+1:N_1+N_2
    b(i,1) = 0.5*c_h(i-N_1,1)*norm(Q_inf)*(Cl_0+Cl_alpha*(AoA_t)+Cl_delta*delta_t);
end
gamma = a\b;
end