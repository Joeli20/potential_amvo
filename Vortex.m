function [V_f_modul,V_x,V_z,Cp_f,CL,Cl,L, CM0, Cm_0, a, gamma] = Vortex(Q_inf, AoA, cosinus,sinus,l_p,nodes,center,normal, tangent)

a = zeros(size(normal,2),size(normal,2));
b = zeros(size(normal,2),1);
gamma = zeros(size(normal,2),1);
Cm_0 = zeros(size(normal,2),1);
trigger1 = size(normal,2);
V_inf(1,1) = Q_inf*cosd(AoA);
V_inf(1,2) = Q_inf*sind(AoA);
CL=0;

for i = 1:size(normal,2)
    for j = 1:size(normal,2)
        b(i,1) = -(V_inf(1,1)*tangent(1,i)+V_inf(1,2)*tangent(2,i));
        if i == j
            a(i,j) = -1/2;
        else
            x_pan(i,j) = (center(1,i)-nodes(1,j))*cosinus(1,j)-(center(2,i)-nodes(2,j))*sinus(1,j);
            z_pan(i,j) = (center(1,i)-nodes(1,j))*sinus(1,j)+(center(2,i)-nodes(2,j))*cosinus(1,j);
            r1(i,j) = sqrt(x_pan(i,j).^2+z_pan(i,j).^2);
            r2(i,j) = sqrt((x_pan(i,j)-l_p(j,1)).^2+z_pan(i,j).^2);
            w_ind_p(i,j) = 1/(4*pi)*log((r2(i,j).^2)/(r1(i,j).^2)); 
            theta1(i,j) = atan2(z_pan(i,j),x_pan(i,j));
            theta2(i,j) = atan2(z_pan(i,j),(x_pan(i,j)-l_p(j,1)));            
            u_ind_p(i,j) = (theta2(i,j)-theta1(i,j))/(2*pi);
            u_ind_global(i,j) = u_ind_p(i,j)*cosinus(1,j) + w_ind_p(i,j)*sinus(1,j);
            w_ind_global(i,j) = (-u_ind_p(i,j))*sinus(1,j) + w_ind_p(i,j)*cosinus(1,j);
            V(i,j) = sqrt(u_ind_global(i,j).^2 +w_ind_global(i,j).^2);
            a(i,j) = u_ind_global(i,j)*cosinus(1,i)+w_ind_global(i,j)*-sinus(1,i); 
        end

    end

end
    k=round(trigger1/4);
    a(k,:)=0;
    a(k,1)=1;
    a(k,trigger1)=1;
    b(k,1) = 0;
    %a(k,k)=-1/2;


    gamma = a\b;
    gamma(k,1) = (gamma(k-1,1)+gamma(k+1,1))/2;

for i = 1:size(normal,2)
    V_x(i) = gamma(i)*tangent(1,i);
    V_z(i) = gamma(i)*tangent(2,i);
    V_f_modul(i) = abs(gamma(i));
    Cp_f(i) = 1-(gamma(i)/Q_inf)^2;
    Cl(i) = (2*gamma(i).*l_p(i,1))./(Q_inf*1);
    L_rara(i) = (gamma(i).*l_p(i,1));
    Cm_0(i) = Cm_0(i)+(Cp_f(i)/1^2).*(((center(1,i)-1/4)*(nodes(1,i+1)-nodes(1,i)))+(center(2,i)*(nodes(2,i+1)-nodes(2,i))));
end
L = sum(L_rara)*Q_inf;
CL = sum(Cl);
CM0 = sum(Cm_0);
end