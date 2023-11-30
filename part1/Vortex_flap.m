function [v_f,v_x,v_z,cp,cl,cm_0,gamma, a] = Vortex_flap(Q_inf,AoA,cosinus_m,cosinus_f,sinus_m,sinus_f,l_p_m,l_p_f,node_m,node_f,center_m,center_f,vec_tm,vec_tf,c)
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
a = zeros(size(vec_tm,2)+size(vec_tf,2),size(vec_tm,2)+size(vec_tf,2));
b = zeros(size(vec_tm,2)+size(vec_tf,2),1);
cm_0_nod = zeros(size(vec_tm,2)+size(vec_tf,2),1);
trigger1 = size(vec_tm,2)+size(vec_tf,2);
trigger2 = size(vec_tm,2);
trigger3 = size(vec_tf,2);
V_inf(1,1) = Q_inf*cosd(AoA);
V_inf(1,2) = Q_inf*sind(AoA);
x_pan = zeros(size(vec_tm,2)+size(vec_tf,2),size(vec_tm,2)+size(vec_tf,2));
z_pan = zeros(size(vec_tm,2)+size(vec_tf,2),size(vec_tm,2)+size(vec_tf,2));
r1 = zeros(size(vec_tm,2)+size(vec_tf,2),size(vec_tm,2)+size(vec_tf,2));
r2 = zeros(size(vec_tm,2)+size(vec_tf,2),size(vec_tm,2)+size(vec_tf,2));
w_ind_p = zeros(size(vec_tm,2)+size(vec_tf,2),size(vec_tm,2)+size(vec_tf,2));
theta1 = zeros(size(vec_tm,2)+size(vec_tf,2),size(vec_tm,2)+size(vec_tf,2));
theta2 = zeros(size(vec_tm,2)+size(vec_tf,2),size(vec_tm,2)+size(vec_tf,2));
u_ind_p = zeros(size(vec_tm,2)+size(vec_tf,2),size(vec_tm,2)+size(vec_tf,2));
u_ind_global = zeros(size(vec_tm,2)+size(vec_tf,2),size(vec_tm,2)+size(vec_tf,2));
w_ind_global = zeros(size(vec_tm,2)+size(vec_tf,2),size(vec_tm,2)+size(vec_tf,2));

for i = 1:size(vec_tm,2)
    for j = 1:size(vec_tm,2)
        b(i,1) = -(V_inf(1,1)*vec_tm(1,i)+V_inf(1,2)*vec_tm(2,i));
        if i == j
            a(i,j) = -1/2;
        else
            x_pan(i,j) = (center_m(1,i)-node_m(1,j))*cosinus_m(1,j)-(center_m(2,i)-node_m(2,j))*sinus_m(1,j);
            z_pan(i,j) = (center_m(1,i)-node_m(1,j))*sinus_m(1,j)+(center_m(2,i)-node_m(2,j))*cosinus_m(1,j);
            r1(i,j) = sqrt(x_pan(i,j).^2+z_pan(i,j).^2);
            r2(i,j) = sqrt((x_pan(i,j)-l_p_m(j,1)).^2+z_pan(i,j).^2);
            w_ind_p(i,j) = 1/(4*pi)*log((r2(i,j).^2)/(r1(i,j).^2)); 
            theta1(i,j) = atan2(z_pan(i,j),x_pan(i,j));
            theta2(i,j) = atan2(z_pan(i,j),(x_pan(i,j)-l_p_m(j,1)));            
            u_ind_p(i,j) = (theta2(i,j)-theta1(i,j))/(2*pi);
            u_ind_global(i,j) = u_ind_p(i,j)*cosinus_m(1,j) + w_ind_p(i,j)*sinus_m(1,j);
            w_ind_global(i,j) = (-u_ind_p(i,j))*sinus_m(1,j) + w_ind_p(i,j)*cosinus_m(1,j);
            a(i,j) = u_ind_global(i,j)*cosinus_m(1,i)+w_ind_global(i,j)*-sinus_m(1,i); 
        end

    end
    for j = size(vec_tm,2)+1:size(vec_tf,2)+size(vec_tm,2)
    %for j = 1:size(vec_tf,2)
        if i == j
            a(i,j) = -1/2;
        else
            x_pan(i,j) = (center_m(1,i)-node_f(1,j-size(vec_tm,2))).*cosinus_f(1,j-size(vec_tm,2))-(center_m(2,i)-node_f(2,j-size(vec_tm,2))).*sinus_f(1,j-size(vec_tm,2));
            z_pan(i,j) = (center_m(1,i)-node_f(1,j-size(vec_tm,2))).*sinus_f(1,j-size(vec_tm,2))+(center_m(2,i)-node_f(2,j-size(vec_tm,2))).*cosinus_f(1,j-size(vec_tm,2));
            r1(i,j) = sqrt(x_pan(i,j).^2+z_pan(i,j).^2);
            r2(i,j) = sqrt((x_pan(i,j)-l_p_f(j-size(vec_tm,2),1)).^2+z_pan(i,j).^2);
            w_ind_p(i,j) = 1/(4*pi)*log((r2(i,j).^2)/(r1(i,j).^2)); 
            theta1(i,j) = atan2(z_pan(i,j),x_pan(i,j));
            theta2(i,j) = atan2(z_pan(i,j),(x_pan(i,j)-l_p_f(j-size(vec_tm,2),1)));            
            u_ind_p(i,j) = (theta2(i,j)-theta1(i,j))/(2*pi);
            u_ind_global(i,j) = u_ind_p(i,j)*cosinus_f(1,j-size(vec_tm,2)) + w_ind_p(i,j)*sinus_f(1,j-size(vec_tm,2));
            w_ind_global(i,j) = (-u_ind_p(i,j))*sinus_f(1,j-size(vec_tm,2)) + w_ind_p(i,j)*cosinus_f(1,j-size(vec_tm,2));
            a(i,j) = u_ind_global(i,j)*cosinus_m(1,i)+w_ind_global(i,j)*-sinus_m(1,i); 
        end

    end

end
for i = size(vec_tm,2)+1:size(vec_tf,2)+size(vec_tm,2)
%for i = 1:size(vec_tf,2)
    b(i,1) = -(V_inf(1,1)*vec_tf(1,i-size(vec_tm,2))+V_inf(1,2)*vec_tf(2,i-size(vec_tm,2)));
    for j = 1:size(vec_tm,2)
        if i == j
            a(i,j) = -1/2;
        else
            x_pan(i,j) = (center_f(1,i-size(vec_tm,2))-node_m(1,j))*cosinus_m(1,j)-(center_f(2,i-size(vec_tm,2))-node_m(2,j))*sinus_m(1,j);
            z_pan(i,j) = (center_f(1,i-size(vec_tm,2))-node_m(1,j))*sinus_m(1,j)+(center_f(2,i-size(vec_tm,2))-node_m(2,j))*cosinus_m(1,j);
            r1(i,j) = sqrt(x_pan(i,j).^2+z_pan(i,j).^2);
            r2(i,j) = sqrt((x_pan(i,j)-l_p_m(j,1)).^2+z_pan(i,j).^2);
            w_ind_p(i,j) = 1/(4*pi)*log((r2(i,j).^2)/(r1(i,j).^2)); 
            theta1(i,j) = atan2(z_pan(i,j),x_pan(i,j));
            theta2(i,j) = atan2(z_pan(i,j),(x_pan(i,j)-l_p_m(j,1)));            
            u_ind_p(i,j) = (theta2(i,j)-theta1(i,j))/(2*pi);
            u_ind_global(i,j) = u_ind_p(i,j)*cosinus_m(1,j) + w_ind_p(i,j)*sinus_m(1,j);
            w_ind_global(i,j) = (-u_ind_p(i,j))*sinus_m(1,j) + w_ind_p(i,j)*cosinus_m(1,j);
            a(i,j) = u_ind_global(i,j)*cosinus_f(1,i-size(vec_tm,2))+w_ind_global(i,j)*-sinus_f(1,i-size(vec_tm,2)); 
        end

    end
    for j = size(vec_tm,2)+1:size(vec_tf,2)+size(vec_tm,2)
        if i == j
            a(i,j) = -1/2;
        else
            x_pan(i,j) = (center_f(1,i-size(vec_tm,2))-node_f(1,j-size(vec_tm,2)))*cosinus_f(1,j-size(vec_tm,2))-(center_f(2,i-size(vec_tm,2))-node_f(2,j-size(vec_tm,2)))*sinus_f(1,j-size(vec_tm,2));
            z_pan(i,j) = (center_f(1,i-size(vec_tm,2))-node_f(1,j-size(vec_tm,2)))*sinus_f(1,j-size(vec_tm,2))+(center_f(2,i-size(vec_tm,2))-node_f(2,j-size(vec_tm,2)))*cosinus_f(1,j-size(vec_tm,2));
            r1(i,j) = sqrt(x_pan(i,j).^2+z_pan(i,j).^2);
            r2(i,j) = sqrt((x_pan(i,j)-l_p_f(j-size(vec_tm,2),1)).^2+z_pan(i,j).^2);
            w_ind_p(i,j) = 1/(4*pi)*log((r2(i,j).^2)/(r1(i,j).^2)); 
            theta1(i,j) = atan2(z_pan(i,j),x_pan(i,j));
            theta2(i,j) = atan2(z_pan(i,j),(x_pan(i,j)-l_p_f(j-size(vec_tm,2),1)));            
            u_ind_p(i,j) = (theta2(i,j)-theta1(i,j))/(2*pi);
            u_ind_global(i,j) = u_ind_p(i,j)*cosinus_f(1,j-size(vec_tm,2)) + w_ind_p(i,j)*sinus_f(1,j-size(vec_tm,2));
            w_ind_global(i,j) = (-u_ind_p(i,j))*sinus_f(1,j-size(vec_tm,2)) + w_ind_p(i,j)*cosinus_f(1,j-size(vec_tm,2));
            a(i,j) = u_ind_global(i,j)*cosinus_f(1,i-size(vec_tm,2))+w_ind_global(i,j)*-sinus_f(1,i-size(vec_tm,2)); 
        end
    end

end
    k=round(trigger2/4);
    a(k,:)=0;
    a(k,1)=1;
    a(k,trigger2)=1;
    b(k,1) = 0;
    g=round((trigger2+(trigger3/4)));
    a(g,:)=0;
    a(g,trigger2+1)=1;
    a(g,trigger1)=1;
    b(g,1) = 0;   
    %a(k,k)=-1/2;


    gamma = a\b;
    gamma(k,1) = (gamma(k-1,1)+gamma(k+1,1))/2;
    gamma(g,1) = (gamma(g-1,1)+gamma(g+1,1))/2;
 
for i = 1:size(vec_tm,2)
    v_x(i) = gamma(i)*vec_tm(1,i);
    v_z(i) = gamma(i)*vec_tm(2,i);
    v_f(i) = abs(gamma(i));
    cp(i) = 1-(gamma(i)/Q_inf)^2;
    cl_nod(i) = (2*gamma(i).*l_p_m(i,1))./(Q_inf*c);
    L_rara(i) = (gamma(i).*l_p_m(i,1));
    cm_0_nod(i) = cm_0_nod(i)+(cp(i)/c^2).*(((center_m(1,i)-c/4)*(node_m(1,i+1)-node_m(1,i)))+(center_m(2,i)*(node_m(2,i+1)-node_m(2,i))));
end
for i = size(vec_tm,2)+1:size(vec_tf,2)+size(vec_tm,2)
    v_x(i) = gamma(i)*vec_tf(1,i-size(vec_tm,2));
    v_z(i) = gamma(i)*vec_tf(2,i-size(vec_tm,2));
    v_f(i) = abs(gamma(i));
    cp(i) = 1-(gamma(i)/Q_inf)^2;
    cl_nod(i) = (2*gamma(i).*l_p_f(i-size(vec_tm,2),1))./(Q_inf*c);
    L_rara(i) = (gamma(i).*l_p_f(i-size(vec_tm,2),1));
    cm_0_nod(i) = cm_0_nod(i)+(cp(i)/c^2).*(((center_f(1,i-size(vec_tm,2))-c/4)*(node_f(1,i+1-size(vec_tm,2))-node_f(1,i-size(vec_tm,2))))+(center_f(2,i-size(vec_tm,2))*(node_f(2,i+1-size(vec_tm,2))-node_f(2,i-size(vec_tm,2)))));
end
L = sum(L_rara)*Q_inf;
cl = sum(cl_nod);
cm_0 = sum(cm_0_nod);
end