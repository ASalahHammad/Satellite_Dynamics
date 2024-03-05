clear
clc
close all

%% Satellite Universat
v = 7.5465e3; % Linear Velocity of Satellite (m/s)
h = 619e3; % Height above Earth (m.)
RE = 6367e3; % Radius of Earth (m.)
w0 = v/(RE+h); % Angular Velocity of Satellite around Earth (rad/s)
Ix = 6.916894; % Kg.m^2
Iy = 6.684942; % Kg.m^2
Iz = 4.915737; % Kg.m^2

sx = (Iy-Iz)/Ix;
sy = (Ix-Iz)/Iy;
sz = abs(Iy-Ix)/Iz;

Dx = 1.0; % Damping Coefficient
Dy = 2.0; % Damping Coefficient
Dz = 1.0; % Damping Coefficient
Iwx = 2; % Kg.m^2
Iwy = 2; % Kg.m^2
Iwz = 2; % Kg.m^2




%% ==================== Simulink ====================
%% ========== Approach 1 ==========
% c5 = ((Dx*Iwx*Iwz*Iz + Dz*Iwx*Iwz*Ix + Dx*Iwz*Ix*Iz + Dz*Iwx*Ix*Iz)/(Iwx*Iwz*Ix*Iz));
% c4 = ((Dx*Dz*Iwx*Iwz + Dx*Dz*Iwz*Ix + Dx*Dz*Iwx*Iz + Dx*Dz*Ix*Iz + Iwx*Iwz*Ix*Iz*w0 - Iwx*Iwz*Ix*Iz*sx*w0 - Iwx*Iwz*Ix*Iz*sz*w0 + 4*Iwx*Iwz*Ix*Iz*sx*w0^2 + Iwx*Iwz*Ix*Iz*sz*w0^2 + Iwx*Iwz*Ix*Iz*sx*sz*w0)/(Iwx*Iwz*Ix*Iz));
% c3 = ((Dx*Iwz*Ix*Iz*w0 + Dz*Iwx*Ix*Iz*w0 - Dx*Iwz*Ix*Iz*sx*w0 - Dx*Iwz*Ix*Iz*sz*w0 - Dz*Iwx*Ix*Iz*sx*w0 - Dz*Iwx*Ix*Iz*sz*w0 + 4*Dz*Iwx*Iwz*Ix*sx*w0^2 + Dx*Iwx*Iwz*Iz*sz*w0^2 + 4*Dx*Iwz*Ix*Iz*sx*w0^2 + Dx*Iwz*Ix*Iz*sz*w0^2 + 4*Dz*Iwx*Ix*Iz*sx*w0^2 + Dz*Iwx*Ix*Iz*sz*w0^2 + Dx*Iwz*Ix*Iz*sx*sz*w0 + Dz*Iwx*Ix*Iz*sx*sz*w0)/(Iwx*Iwz*Ix*Iz));
% c2 = ((Dx*Dz*Ix*Iz*w0 - Dx*Dz*Ix*Iz*sx*w0 - Dx*Dz*Ix*Iz*sz*w0 + 4*Dx*Dz*Iwz*Ix*sx*w0^2 + Dx*Dz*Iwx*Iz*sz*w0^2 + 4*Dx*Dz*Ix*Iz*sx*w0^2 + Dx*Dz*Ix*Iz*sz*w0^2 + 4*Iwx*Iwz*Ix*Iz*sx*sz*w0^4 + Dx*Dz*Ix*Iz*sx*sz*w0)/(Iwx*Iwz*Ix*Iz));
% c1 = ((4*Dx*Iwz*Ix*Iz*sx*sz*w0^4 + 4*Dz*Iwx*Ix*Iz*sx*sz*w0^4)/(Iwx*Iwz*Ix*Iz));
% c0 = (4*Dx*Dz*sx*sz*w0^4)/(Iwx*Iwz);
% 
% c5 = ((Dx*Iwx*Iwz*Iz^2 + Dx*Iwz*Ix*Iz^2 + Dz*Iwx*Ix*Iz^2 + Dz*Iwx*Iwz*Ix*Iz)/(Iwx*Iwz*Ix*Iz^2));
% c4 = ((Dx*Dz*Iwx*Iz^2 + Dx*Dz*Ix*Iz^2 - Dz^2*Iwx*Iwz*Ix + Dz^2*Iwx*Ix*Iz + Dx*Dz*Iwx*Iwz*Iz + Dx*Dz*Iwz*Ix*Iz + 4*Iwx*Iwz*Ix*Iz^2*sx*w0^2 + Iwx*Iwz*Ix*Iz^2*sz*w0^2)/(Iwx*Iwz*Ix*Iz^2));
% c3 = ((Dx*Dz^2*Iwx*Iz - Dx*Dz^2*Iwz*Ix - Dx*Dz^2*Iwx*Iwz + Dx*Dz^2*Ix*Iz + Iwx*Iwz*Ix*Iz^2*w0^2 + 4*Dx*Iwz*Ix*Iz^2*sx*w0^2 + Dx*Iwx*Iwz*Iz^2*sz*w0^2 + 4*Dz*Iwx*Ix*Iz^2*sx*w0^2 + Dx*Iwz*Ix*Iz^2*sz*w0^2 + Dz*Iwx*Ix*Iz^2*sz*w0^2 - Iwx*Iwz*Ix*Iz^2*sx*w0^2 - Iwx*Iwz*Ix*Iz^2*sz*w0^2 + 4*Dz*Iwx*Iwz*Ix*Iz*sx*w0^2 + Iwx*Iwz*Ix*Iz^2*sx*sz*w0^2)/(Iwx*Iwz*Ix*Iz^2));
% c2 = ((Dx*Iwz*Ix*Iz^2*w0^2 + Dz*Iwx*Ix*Iz^2*w0^2 + 4*Dx*Dz*Ix*Iz^2*sx*w0^2 + Dx*Dz*Iwx*Iz^2*sz*w0^2 + Dx*Dz*Ix*Iz^2*sz*w0^2 - 4*Dz^2*Iwx*Iwz*Ix*sx*w0^2 - Dx*Iwz*Ix*Iz^2*sx*w0^2 - Dz*Iwx*Ix*Iz^2*sx*w0^2 + 4*Dz^2*Iwx*Ix*Iz*sx*w0^2 - Dx*Iwz*Ix*Iz^2*sz*w0^2 - Dz*Iwx*Ix*Iz^2*sz*w0^2 + 4*Dx*Dz*Iwz*Ix*Iz*sx*w0^2 + Dx*Iwz*Ix*Iz^2*sx*sz*w0^2 + Dz*Iwx*Ix*Iz^2*sx*sz*w0^2 + 4*Iwx*Iwz*Ix*Iz^2*sx*sz*w0^4)/(Iwx*Iwz*Ix*Iz^2));
% c1 = ((Dx*Dz*Ix*Iz^2*w0^2 - 4*Dx*Dz^2*Iwz*Ix*sx*w0^2 - Dx*Dz*Ix*Iz^2*sx*w0^2 + 4*Dx*Dz^2*Ix*Iz*sx*w0^2 - Dx*Dz*Ix*Iz^2*sz*w0^2 + Dx*Dz*Ix*Iz^2*sx*sz*w0^2 + 4*Dx*Iwz*Ix*Iz^2*sx*sz*w0^4 + 4*Dz*Iwx*Ix*Iz^2*sx*sz*w0^4)/(Iwx*Iwz*Ix*Iz^2));
% c0 = (4*Dx*Dz*sx*sz*w0^4)/(Iwx*Iwz);
% sim("sat_sim_1.slx")

%% ========== Approach 2 not complete ==========
% A = [0               1             0               0               0               0               0               0               0;
%      -4*w0^2*sx   -Dx/Ix         Dx/Ix             0               0               0           w0*(1-sx)           0               0;
%      0             Dx/Iwx       -Dx/Iwx            0               0               0               0               0               0;
%      0               0             0               0               1               0               0               0               0;
%      0               0             0          -3*w0^2*sy        -Dy/Iy           Dy/Iy             0               0               0;
%      0               0             0               0             Dy/Iwy         -Dy/Iwy            0               0               0;
%      0               0             0               0               0               0               0               1               0;
%      0            -w0*(1-sz)       0               0               0               0           -w0^2*sz         -Dz/Iz         Dz/Iz;
%      0               0             0               0               0               0               0             Dz/Iwz      -Dz/Iwz;
%      ];
% B = [0 0 0; 1/Ix 0 0; 0 0 0; 0 0 0; 0 1/Iy 0; 0 0 0; 0 0 0; 0 0 1/Iz; 0 0 0];
% [NUM, DEN] = ss2tf(A, B, A, B, 1); %% last input iu indicates the iu'th input signal
% TF = cell(1, 9);
% num = cell(9, 10);
% for i=1:9
%     num(i) = {NUM(i, :)};
%     TF(i) = tf(num1, DEN);
% end
% 

%% ========== State Space ========== 
A = [0               1             0               0               0               0               0               0               0;
     -4*w0^2*sx   -Dx/Ix         Dx/Ix             0               0               0            w0*(1-sx)           0               0;
     0             Dx/Iwx       -Dx/Iwx            0               0               0               0               0               0;
     0               0             0               0               1               0               0               0               0;
     0               0             0          -3*w0^2*sy        -Dy/Iy           Dy/Iy             0               0               0;
     0               0             0               0             Dy/Iwy         -Dy/Iwy            0               0               0;
     0               0             0               0               0               0               0               1               0;
     0            -w0*(1-sz)       0               0               0               0           -w0^2*sz         -Dz/Iz         Dz/Iz;
     0               0             0               0               0               0               0             Dz/Iwz      -Dz/Iwz;
     ];
B = [0 0 0; 1/Ix 0 0; 0 0 0; 0 0 0; 0 1/Iy 0; 0 0 0; 0 0 0; 0 0 1/Iz; 0 0 0];
C = eye(9);
D = zeros(9, 3);
x0 = [0 0 0 0 0 0 0 0 0].'; % initial conditions
t = 0:100:5e7;
u = 1 * ones(3,numel(t)); % input
% u = [zeros(1, 100) 1000, 0*ones(1, numel(t)-100-1)] .* [1; 1; 1];
[Y, X] = lsim(A, B, C, D, u, t, x0);
plot(t, Y(:, 4)); % theta
hold on; grid on;
% plot(t, Y(:, 1)); % phi
% plot(t, Y(:, 7)); % psi
legend("\theta", "\phi", "\psi")
% sys = ss2tf(A, B, A, B, 1)


