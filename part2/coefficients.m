function [L_total, L_tail,Cl_c_wing,L_wing,CL_wing, CL_tail, Cl_c_tail, CD_v_wing, Cd_wing, Cd_tail, CD_v_tail,CD_i_wing, CD_i_tail, CD_wing,CD_tail, Cd_i_wing,Cd_i_tail, Cd_v_wing,Cd_v_tail, alpha_i_wing,alpha_i_tail, Cm_0_ala,Cm_0_tail,Effi_ala,Effi_tail,Cm_roll_ala,Cm_roll_tail, M_roll_ala,M_roll_tail,N_mom_cua,C_n_cua] = coefficients(rho,y,y_h,x_c,y_c,z_c,x_c_h,y_c_h,z_c_h,c,c_h,Q_inf,AoA, AoA_t,N_1,N_2,Cl_0,Cl_alpha,gamma, b_wing,b_tail, S_wing, S_tail,theta,Cd_v0,Cd_v_Cl,M_chord_ala,M_chord_tail, Cm_025)
%Ratio_wing
% Cm_0_wing, Cm_0_tail, CM_roll_wing, CM_roll_tail

for i = 1:N_1
    L(i,1) = rho*norm(Q_inf)*gamma(i,1)*(y(1,i+1)-y(1,i));
    CL(i,1) = (2/(norm(Q_inf)*S_wing))*gamma(i,1)*(y(1,i+1)-y(1,i));
    Cl_c_wing(i,1) = (2*gamma(i,1))/(c(i,1)*norm(Q_inf));
    alpha_i_wing(i,1) = (((Cl_c_wing(i,1)-Cl_0)/Cl_alpha)-AoA-theta(i,1));
    D_i_wing(i,1) = -rho*norm(Q_inf)*gamma(i,1)*(y(1,i+1)-y(1,i))*alpha_i_wing(i,1);
    Cd_i_wing(i,1) = -2*((gamma(i,1)*(y(1,i+1)-y(1,i))*alpha_i_wing(i,1))/(norm(Q_inf)*S_wing));
    Cd_v_wing(i,1) = Cd_v0+Cd_v_Cl*(Cl_c_wing(i,1))^2;
    Cd_v_wing2(i,1) = Cd_v_wing(i,1)*c(i,1)*(y(1,i+1)-y(1,i));
    D_v_wing(i,1) = 0.5*rho*(norm(Q_inf))^2*Cd_v_wing2(i,1);
    Cd_wing(i,1) = Cd_i_wing(i,1)+Cd_v_wing(i,1);
    Cm_x(i,1) = -2*(x_c(i,1)*gamma(i,1)*(y(1,i+1)-y(1,i)))/(norm(Q_inf)*S_wing*M_chord_ala);
    Cm_y(i,1) = -2*(y_c(i,1)*gamma(i,1)*(y(1,i+1)-y(1,i)))/(norm(Q_inf)*S_wing*M_chord_ala);
    Cm_z(i,1) = -2*(z_c(i,1)*gamma(i,1)*(y(1,i+1)-y(1,i)))/(norm(Q_inf)*S_wing*M_chord_ala);
    Cm_y2(i,1) = (y_c(i,1)*gamma(i,1)*(y(1,i+1)-y(1,i)));
end

L_wing = sum(L);
CL_wing = sum(CL);
CD_i_wing = sum(Cd_i_wing);
CD_v_wing = sum(Cd_v_wing)/S_wing;
CD_wing = CD_i_wing+CD_v_wing;
D_v_w = sum(D_v_wing);
D_i_w = sum(D_i_wing);
D_total_w = D_v_w+D_i_w;
CM_x = sum(Cm_x);
CM_y = sum(Cm_y);
CM_z = sum(Cm_z);
Cm_0_ala = Cm_025 - CM_x - CM_y - CM_z;
Effi_ala = CL_wing/CD_wing;
M_roll_ala = sum(Cm_y2)*rho*norm(Q_inf);
Cm_roll_ala = M_roll_ala/(0.5*rho*(norm(Q_inf))^2*S_wing*b_wing);

for i = N_1+1:N_1+N_2
    L(i,1) = rho*norm(Q_inf)*gamma(i,1)*(y_h(1,i+1-N_1)-y_h(1,i-N_1));
    l_tail(i-N_1,1) = rho*norm(Q_inf)*gamma(i,1)*(y_h(1,i+1-N_1)-y_h(1,i-N_1));
    CL(i,1) = (2/(norm(Q_inf)*S_tail))*gamma(i,1)*(y_h(1,i+1-N_1)-y_h(1,i-N_1));
    CL_tail(i-N_1,1) = (2/(norm(Q_inf)*S_tail))*gamma(i,1)*(y_h(1,i+1-N_1)-y_h(1,i-N_1));
    Cl_c_tail(i-N_1,1) = (2*gamma(i,1))/(c_h(i-N_1,1)*norm(Q_inf));
    alpha_i_tail(i-N_1,1) = ((Cl_c_tail(i-N_1,1)-Cl_0)/Cl_alpha)-AoA_t;
    D_i_tail(i-N_1,1) = -rho*norm(Q_inf)*gamma(i,1)*(y_h(1,i+1-N_1)-y_h(1,i-N_1))*alpha_i_tail(i-N_1,1);
    Cd_i_tail(i-N_1,1) = -2*((gamma(i,1)*(y_h(1,i+1-N_1)-y_h(1,i-N_1))*alpha_i_tail(i-N_1,1))/(norm(Q_inf)*S_tail));
    Cd_v_tail(i-N_1,1) = Cd_v0+Cd_v_Cl*(Cl_c_tail(i-N_1,1))^2;
    Cd_v_tail2(i-N_1,1) = Cd_v_tail(i-N_1,1)*c_h(i-N_1,1)*(y_h(1,i+1-N_1)-y_h(1,i-N_1));
    D_v_tail(i-N_1,1) = 0.5*rho*(norm(Q_inf))^2*Cd_v_tail2(i-N_1,1);
    Cd_tail(i-N_1,1) = Cd_i_tail(i-N_1,1)+Cd_v_tail(i-N_1,1);
    Cm_x_h(i-N_1,1) = -2*(x_c_h(i-N_1,1)*gamma(i,1)*(y(1,i+1-N_1)-y(1,i-N_1)))/(norm(Q_inf)*S_tail*M_chord_tail);
    Cm_y_h(i-N_1,1) = -2*(y_c_h(i-N_1,1)*gamma(i,1)*(y(1,i+1-N_1)-y(1,i-N_1)))/(norm(Q_inf)*S_tail*M_chord_tail);
    Cm_z_h(i-N_1,1) = -2*(z_c_h(i-N_1,1)*gamma(i,1)*(y(1,i+1-N_1)-y(1,i-N_1)))/(norm(Q_inf)*S_tail*M_chord_tail);

end
L_tail = sum(l_tail);
CL_tail = sum(CL_tail);
L_total = L_tail+L_wing;
CD_i_tail = sum(Cd_i_tail);
CD_v_tail = sum(Cd_v_tail2)/S_tail;
CD_tail = CD_i_tail+CD_v_tail;
D_v_t = sum(D_v_tail);
D_i_t = sum(D_i_tail);
D_total_t = D_v_t+D_i_t;
CM_x_h = sum(Cm_x_h);
CM_y_h = sum(Cm_y_h);
CM_z_h = sum(Cm_z_h);
Cm_0_tail = Cm_025 - CM_x_h - CM_y_h - CM_z_h;
Effi_tail = CL_tail/CD_tail;
M_roll_tail = rho*norm(Q_inf)*sum(Cm_y_h);
Cm_roll_tail = M_roll_tail/(0.5*rho*(norm(Q_inf))^2*S_tail*b_tail);
for i = 1:N_2
    N_mom_cua(i,1) = -(y_c_h(i,1)*(D_v_tail(i,1)+D_i_tail(i,1)));
end
C_n_cua = -sum(N_mom_cua)/(0.5*rho*(norm(Q_inf))^2*S_tail*b_tail);
end