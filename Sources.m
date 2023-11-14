function [v_f,cp,a,sigma] = Sources(Q_inf,AoA,cosinus,sinus,l_p,node,center,vec_n)
% Sources strenght calculation
%
% Written by: Joel Campo, Jordi Gallart, Martí Santamaria, 2023
% Group 16. AMVO. MUEA.
%
% Inputs:
%   Q_inf: Fluid velocity [m/s]
%   AoA: Angle of attack [degrees]
%   cosinus: Cosinus of every panel
%   sinus: Sinus of every panel
%   l_p: Lenght of every panel
%   node: Position of every node
%   center: Position of the center point of a panel
%   vec_n: Normal vector
% Outputs:
%   v_f: Module of final velocity
%   cp: Pressure coefficient
%   sigma: Sources strenght
%   a: DELETE

% Preallocating
a = zeros(size(vec_n,2),size(vec_n,2));
b = zeros(size(vec_n,2),1);
valor_sum=zeros(size(vec_n,2),1);

for i = 1:size(vec_n,2)
    for j = 1:size(vec_n,2)
        b(i,1) = -Q_inf*vec_n(1,i);
        if i == j
            a(i,j) = 1/2;
        else
            x_pan(i,j) = (center(1,i)-node(1,j))*cosinus(1,j)-(center(2,i)-node(2,j))*sinus(1,j);
            z_pan(i,j) = (center(1,i)-node(1,j))*sinus(1,j)+(center(2,i)-node(2,j))*cosinus(1,j);
            r1(i,j) = sqrt(x_pan(i,j).^2+z_pan(i,j).^2);
            r2(i,j) = sqrt((x_pan(i,j)-l_p(j,1)).^2+z_pan(i,j).^2);
            u_ind_p(i,j) = 1/(4*pi)*log((r1(i,j).^2)/(r2(i,j).^2)); 
            theta1(i,j) = atan2(z_pan(i,j),x_pan(i,j));
            theta2(i,j) = atan2(z_pan(i,j),(x_pan(i,j)-l_p(j,1)));
            w_ind_p(i,j) = (theta2(i,j)-theta1(i,j))/(2*pi);
            u_ind_global(i,j) = u_ind_p(i,j)*cosinus(1,j) + w_ind_p(i,j)*sinus(1,j);
            w_ind_global(i,j) = -u_ind_p(i,j)*sinus(1,j) + w_ind_p(i,j)*cosinus(1,j);
            V(i,j) = sqrt(u_ind_global(i,j).^2 +w_ind_global(i,j).^2);
            a(i,j) = u_ind_global(i,j)*sinus(1,i)+w_ind_global(i,j)*cosinus(1,i);
        end

    end

end
    sigma = a\b;

for i = 1:size(vec_n,2)
    for j = 1:size(vec_n,2)
        valor_sum(i) = valor_sum(i)+sigma(i)*V(i,j);
    end
    v_f(i) = Q_inf+valor_sum(i);
    cp(i) = 1-((v_f(i)).^2)./(Q_inf^2);
end
end