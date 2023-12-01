clc;
clear;
close all;

%% Input values
% Wing Geometry
b = 4.5; %[m] wingspan
c_r = 0.9; %[m] chord at the root
c_t = 0.6; %[m] chord at the tip
AoA = 0; %[º] Angle of attack

% Other Wing values
m = 40; %[kg] mass of the wing
rho = 1.225; %[kg/m^3] air density
g = 9.81; %[m/s^2] acceleration of gravity

Cl_0 = 0.24;
Cl_a = 6.7;

%% Previous calculations
S = b * ((c_r + c_t)/2);
W = m*g;

Cl = Cl_0 + Cl_a*AoA;

%% Main matter

U_inf = sqrt((2*W)/(rho*S*CL));