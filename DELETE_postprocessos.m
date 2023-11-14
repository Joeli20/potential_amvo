clear;
clc;
close all;

%% 1
% attack_angles = [2,4,6,8,10];
% for i = 1:length(attack_angles)
    % "cosetes"
    % rho, Q_inf, c, vortex_strength(i), panel_length(i)
    
    % LIFT COEFFICIENT
    % Cl = zeros(1,length(attack_angles));
    % sumat(i) = sum((vortex_strength .* panel_length)/(Q_inf*c));
    % Cl(i) = 2*sumat(i);

    % PITCHING MOMENT COEFFICIENT 1/4
    % Cp = 1-(vortex_strength/Q_inf).^2;
    % Cm1_4(i) = sum((Cp/c^2).*((x_c-c/4).*(diff(x))+z_c.*(diff(z)));

% end

%% 2
% attack_angles = [0,2,4,6];
% for i = 1:length(attack_angles)
    % "cosetes"
    % Mcr = sym('Mcr');
    % Cp = 1-(vortex_strength/Q_inf).^2;
    % Cp0 = min(Cp);
    % Cp = Cp0/(sqrt(1-Mcr^2)+(Cp0/2)*(Mcr^2/sqrt(1-Mcr^2))*(1+(gamma-1)*Mcr^2/2));
    % Cp_cr = (2/(gamma*Mcr^2))*(((2+(gamma-1)*Mcr^2)/(1+gamma))^(gamma/(gamma-1))-1);
    % Mcritic_eq = Cp-Cp_cr == 0;
    % Mcritic_sol(i) = solve(Mcritic_eq, Mcr);

% end

%% 3 (attack_angle = 2)
% Mcr2 = Mcritic_sol(2);
% Minf1 = Mcr2-0.15;
% Minf2 = Mcr2-0.10;
% Minf3 = Mcr2-0.05;
% Minf4 = Mcr2;
% Minf = [Minf1,Minf2,Minf3,Minf4];
% Cl0 = Cl(1);
% for i = 1:length(Minf)
    % beta(i) = sqrt(1-Minf(i)^2);
    % Cl(i) = 1/beta(i) *Cl0;

% end


%% 4
% attack_angle = 0;
% c = c1+c2+d;
% deflect_anlges = [4,8,12,16,20];
% for i = 1:length(deflect_angles)
    % "cosetes"
    % rho, Q_inf, c, vortex_strength(i), panel_length(i)

    % LIFT COEFFICIENT
    % Cl = zeros(1,length(deflect_angles));
    % sumat(i) = sum((vortex_strength .* panel_length)/(Q_inf*c));
    % Cl(i) = 2*sumat(i);

    % PITCHING MOMENT COEFFICIENT 1/4
    % Cp = 1-(vortex_strength/Q_inf).^2;
    % Cm1_4(i) = sum((Cp/c^2).*((x_c-c/4).*(diff(x))+z_c.*(diff(z)));

% end