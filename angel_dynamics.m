function x_dot = angel_dynamics(t, x, u_c, sys, w_const)

%% Primary States
rho = x(1:3); 
nu = x(4:6);
phi = x(7);
theta = x(8);
psi = x(9);    
omega = x(10:12); 
alpha_w = x(13:15);
alpha_w_dot = x(16:18);

%% Control States
T1 = u_c(1);
T2 = u_c(2);
T3 = u_c(3);
W1 = u_c(4);
W2 = u_c(5);
W3 = u_c(6);

%% System Parameters
g = sys(1); 
m = sys(2); 
Jn = sys(3); 
Je = sys(4); 
Jd = sys(5); 
Iw = sys(6); 
k_M = sys(7); 
R_f = sys(8);
J = [Jn,  0,   0; 
      0, Je,   0; 
      0,  0,  Jd]; 
C_NB = C(3, psi) * C(2, theta) * C(1, phi); %DCM to convert from body to inertial
D = D_matrix(phi, theta);
l = w_const(1); 
S = w_const(2); 
y_MAC = w_const(3);
MAC = w_const(4);

%% Arm Rotations and Body Rotation

C_flip = [1 0 0;
          0 1 0;
          0 0 -1];
C_B1 = C(3,pi);      % Arm1 to Body
C_B1([2,3,4]) = 0;
C_B2 = C(3,-pi/3);   % Arm2 to Body
C_B3 = C(3,pi/3);   % Arm3 to Body
C_1P = C(1,alpha_w(1)) * C_flip;   % Prop1 to Arm1
C_2P = C(1,alpha_w(2)) * C_flip;   % Prop2 to Arm2
C_3P = C(1,alpha_w(3)) * C_flip;   % Prop3 to Arm3
C_B1P =  C_B1 * C_1P;       % Prop1 to Body
C_B2P =  C_B2 * C_2P;       % Prop2 to Body
C_B3P =  C_B3 * C_3P;       % Prop3 to Body

%% Arm Length Vectors

l_ = [l; 0; 0];
l1_B = C_B1 * l_;
l2_B = C_B2 * l_;
l3_B = C_B3 * l_;

%% y_MAC Length Vectors

y_ = [y_MAC; 0; 0];
y1_B = C_B1 * y_;
y2_B = C_B2 * y_;
y3_B = C_B3 * y_;

%% MAC Length Vectors

mac_ = [0; 0; -MAC/2]; % Assuming pivot axis is at half the MAC
mac_1P = C_1P * mac_;
mac_2P = C_2P * mac_;
mac_3P = C_3P * mac_;

%% Lift and Drag using formula

% Velocity of each wing at y_MAC due to omega and nu
v1_B = nu + cross(omega, y1_B);
v2_B = nu + cross(omega, y2_B);
v3_B = nu + cross(omega, y3_B);
v1_1 = transpose(C_B1P) * v1_B;
v2_2 = transpose(C_B2P) * v2_B;
v3_3 = transpose(C_B3P) * v3_B;
v1_MAC = norm(v1_1);
v2_MAC = norm(v2_2);
v3_MAC = norm(v3_3);
alpha1_eff = alpha_eff(v1_1(1), v1_1(2),v1_1(3));
alpha2_eff = alpha_eff(v2_2(1), v2_2(2),v2_2(3));
alpha3_eff = alpha_eff(v3_3(1), v3_3(2),v3_3(3));
L1 = (1/2) * CLx(alpha1_eff) * S * v1_MAC^2;
L2 = (1/2) * CLx(alpha2_eff) * S * v2_MAC^2;
L3 = (1/2) * CLx(alpha3_eff) * S * v3_MAC^2;
D1 = (1/2) * CDx(alpha1_eff) * S * v1_MAC^2;
D2 = (1/2) * CDx(alpha2_eff) * S * v2_MAC^2;
D3 = (1/2) * CDx(alpha3_eff) * S * v3_MAC^2;
m1 = (1/2) * Cmx(alpha1_eff) * S * MAC * v1_MAC^2;
m2 = (1/2) * Cmx(alpha2_eff) * S * MAC * v2_MAC^2;
m3 = (1/2) * Cmx(alpha3_eff) * S * MAC * v3_MAC^2;

%% Propeller Thrusts

T1_1P = [0; 0; T1];
T2_2P = [0; 0; T2];
T3_3P = [0; 0; T3];
T1_1 = C_1P * T1_1P;
T2_2 = C_2P * T2_2P;
T3_3 = C_3P * T3_3P;
T1_B = C_B1P * T1_1P;
T2_B = C_B2P * T2_2P;
T3_B = C_B3P * T3_3P;

%% Propeller Motor Moments

M = @(T) k_M * T; % k_M can be negative if counter-clockwise rotation
M1_1P = M(T1_1P);
M2_2P = -M(T2_2P);
M3_3P = -M(T3_3P);
M1_B = C_B1P * M1_1P;
M2_B = C_B2P * M2_2P;
M3_B = C_B3P * M3_3P;

%% Wing Lifts
L1_1P = [0; -L1; 0];
L2_2P = [0; -L2; 0];
L3_3P = [0; -L3; 0];
L1_B = C_B1P * C(1, -alpha1_eff) * L1_1P;
L2_B = C_B2P * C(1, -alpha2_eff) * L2_2P;
L3_B = C_B3P * C(1, -alpha3_eff) * L3_3P;

%% Wing Drags
D1_1P = [0; 0; -D1];
D2_2P = [0; 0; -D2];
D3_3P = [0; 0; -D3];
D1_B = C_B1P * C(1, -alpha1_eff) * D1_1P;
D2_B = C_B2P * C(1, -alpha2_eff) * D2_2P;
D3_B = C_B3P * C(1, -alpha3_eff) * D3_3P;

%% Wing Pitch Moments
m1_1P = [-m1; 0; 0];
m2_2P = [-m2; 0; 0];
m3_3P = [-m3; 0; 0];
m1_B = C_B1P * C(1, -alpha1_eff) * m1_1P;
m2_B = C_B2P * C(1, -alpha2_eff) * m2_2P;
m3_B = C_B3P * C(1, -alpha3_eff) * m3_3P;
m_c = [m1; m2; m3];

%% Wing Blade Motor Moments

W1_1 = [-W1; 0; 0];
W2_2 = [-W2; 0; 0];
W3_3 = [-W3; 0; 0];
W1_B = C_B1 * W1_1;
W2_B = C_B2 * W2_2;
W3_B = C_B3 * W3_3;
W_c = [W1; W2; W3];

%% Resultant Forces in Body Frame

N_f_g = [0; 0; m * g - R_f];      % Total Gravity Force
B_f_T = T1_B + T2_B + T3_B; % Total Thrust Force
B_f_L = L1_B + L2_B + L3_B; % Total Lift Force
B_f_D = D1_B + D2_B + D3_B; % Total Drag Force
B_f_R = transpose(C_NB) * N_f_g + B_f_T + B_f_L + B_f_D;

%% Resultant Torque in Body Frame

B_tau_T = cross(l1_B, T1_B)...
            + cross(l2_B, T2_B) ...
            + cross(l3_B, T3_B);
B_tau_L = cross(y1_B, L1_B)...
            + cross(y2_B, L2_B) ...
            + cross(y3_B, L3_B);
B_tau_D = cross(y1_B, D1_B)...
            + cross(y2_B, D2_B) ...
            + cross(y3_B, D3_B);
B_tau_M = M1_B + M2_B + M3_B;
B_tau_m = m1_B + m2_B + m3_B;
B_tau_W = W1_B + W2_B + W3_B;
B_tau_R = B_tau_T + B_tau_L + B_tau_D + B_tau_M + B_tau_m + B_tau_W;
% B_tau_R = B_tau_T + B_tau_L + B_tau_D + B_tau_M + B_tau_W; % Not considering pitching moments

%% Torque on the Wings Pitch Moment, and Servo Torque
W_m = W_c + m_c; 
% W_m = W_c; 

%% Equations of Motion

% Rate Change States
rho_dot = C_NB * nu;
nu_dot = - cross(omega, nu) + (1/m) * B_f_R;
Lambda_dot = D * omega;
omega_dot =  J \ (B_tau_R - cross(omega, J * omega));
alpha_w_ddot = (1 / Iw) * W_m;

% Dynamics State Function
x_dot = [rho_dot; nu_dot; Lambda_dot; omega_dot; alpha_w_dot; alpha_w_ddot];
