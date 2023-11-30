function [v_f,v_x,v_z,cp, cp_laitone,cl, cl_2,cm_0,gamma] = Vortex2(M_inf,AoA,cosinus,sinus,l_p,node,center,vec_t,vec_n,c)
% Vortex strenght calculation
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
%   vec_t: Tangent vector
%   c: Chord
% Outputs:
%   v_f: Module of final velocity
%   v_x: X component of the Velocity 
%   v_z: Z component of the Velocity 
%   cp: Pressure coefficient
%   cl: Lift coefficient
%   cm_0: Free moment coefficient
%   gamma: Vortex strenght

% Preallocating
a = zeros(size(vec_t,2),size(vec_t,2));
b = zeros(size(vec_t,2),1);
cm_0_nod = zeros(size(vec_t,2),1);
trigger1 = size(vec_t,2);
Q_inf = 343*M_inf;
V_inf(1,1) = Q_inf*cosd(AoA);
V_inf(1,2) = Q_inf*sind(AoA);
GAMMA = 1.4;

for i = 1:size(vec_t,2)
    for j = 1:size(vec_t,2)
        b(i,1) = -(V_inf(1,1)*vec_t(1,i)+V_inf(1,2)*vec_t(2,i));
        if i == j
            a(i,j) = -1/2;
        else
            x_pan(i,j) = (center(1,i)-node(1,j))*cosinus(1,j)-(center(2,i)-node(2,j))*sinus(1,j);
            z_pan(i,j) = (center(1,i)-node(1,j))*sinus(1,j)+(center(2,i)-node(2,j))*cosinus(1,j);
            r1(i,j) = sqrt(x_pan(i,j).^2+z_pan(i,j).^2);
            r2(i,j) = sqrt((x_pan(i,j)-l_p(j,1)).^2+z_pan(i,j).^2);
            w_ind_p(i,j) = 1/(4*pi)*log((r2(i,j).^2)/(r1(i,j).^2)); 
            theta1(i,j) = atan2(z_pan(i,j),x_pan(i,j));
            theta2(i,j) = atan2(z_pan(i,j),(x_pan(i,j)-l_p(j,1)));            
            u_ind_p(i,j) = (theta2(i,j)-theta1(i,j))/(2*pi);
            u_ind_global(i,j) = u_ind_p(i,j)*cosinus(1,j) + w_ind_p(i,j)*sinus(1,j);
            w_ind_global(i,j) = (-u_ind_p(i,j))*sinus(1,j) + w_ind_p(i,j)*cosinus(1,j);
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

for i = 1:size(vec_t,2)
    v_x(i) = gamma(i)*vec_t(1,i);
    v_z(i) = gamma(i)*vec_t(2,i);
    v_f(i) = abs(gamma(i));
    cp(i) = 1-(gamma(i)/Q_inf)^2;
    cp_laitone(i) = cp(i)/(sqrt(1-M_inf^2)+(cp(i)/2)*(M_inf^2/sqrt(1-M_inf^2))*(1+(GAMMA-1)*M_inf^2/2));
    cl_nod(i) = (2*gamma(i).*l_p(i,1))./(Q_inf*c);
    L_rara(i) = (gamma(i).*l_p(i,1));
    cm_0_nod(i) = cm_0_nod(i)+(cp(i)/c^2).*(((center(1,i)-1/4)*(node(1,i+1)-node(1,i)))+(center(2,i)*(node(2,i+1)-node(2,i))));
    cp_laitone2(i) = cp_laitone(i)*l_p(i,1)*vec_n(2,i)/c;
end
cl_22 = sum(cp_laitone2);
cl_2 = -cl_22*cosd(AoA);
L = sum(L_rara)*Q_inf;
cl = sum(cl_nod);
cm_0 = sum(cm_0_nod);
end