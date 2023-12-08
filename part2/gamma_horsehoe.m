function [gamma, V_ij] = gamma_horsehoe(x_c,c,Vinf1,Vinf2,V_vortex_line,N,AoA,Q_inf,Cl_0,Cl_alpha,theta)
V_ij = zeros(N-1,N-1,3);
a = zeros(N-1,N-1);
b = zeros(N-1,1);
for i = 1:size(x_c,1)
    for j = 1:size(x_c,1)
        V_ij(i,j,:) = V_vortex_line(i,j,:)+Vinf1(i,j,:)-Vinf2(i,j,:);
        if i == j
            a(i,j) = -0.5*Cl_alpha*c(i,1)*V_ij(i,j,3)+1;
        else
            a(i,j) = -0.5*Cl_alpha*c(i,1)*V_ij(i,j,3);
        end
    end
    b(i,1) = 0.5*c(i,1)*norm(Q_inf)*(Cl_0+Cl_alpha*(AoA+theta(i,1)));
end
gamma = a\b;
end