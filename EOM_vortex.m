%% State Variables

syms t m Jx Jy Jz g Iw
syms rho_n rho_e rho_d u v w phi theta psi p q r 
syms kp_phi kp_theta kp_psi kd_phi kd_theta kd_psi

%% Symbols for Forces, Moments, Control Inputs

syms alpha1 alpha2 alpha3 alpha % angle of attack with horizontal
syms alpha_dot_1 alpha_dot_2 alpha_dot_3 % rate of change of alpha
% syms Cl_alpha_0 Cl_alpha Cd0 K % Lift and Drag Coefficient elements
syms CL_alpha_0  CL_alpha % Zero lift Intercept and Slope 
syms beta(alpha) % angle of attack with vertical
syms l S y_MAC % arm length to tip motor, Wing area, Wing location of MAC 
% syms Cl1 Cl2 Cl3 Cd1 Cd2 Cd3 % Individual coefficients of Lift and Drag
syms LDR(alpha) % Lift to Drag ratio
syms Cl(alpha) Cd(alpha) % Airfoil coefficient of lift and drag
syms CL(alpha) CD(alpha) % Wing coefficient of lift and drag
syms T1 T2 T3 T % Thrust for each propeller
syms M(T) k_M % Moments for each propeller rotation proportional to Thrust
syms W1 W2 W3 W % Moments for each wing rotation
syms delta1 delta2 delta3 % tip motor inputs
syms sigma1 sigma2 sigma3 % servo motor inputs
syms k1 k2 k3 % Gains mapping to Thrust, Motor Moment, and Wing Moment
syms R_f % Reaction force in vertical direction

%% System Parameters

J_vec= [Jx; Jy; Jz]; 
J = [Jx,  0,   0; 
      0, Jy,   0; 
      0,  0,  Jz]; 
C_NB = C(3, psi) * C(2, theta) * C(1, phi); %DCM to convert from body to inertial
D = D_matrix(phi, theta);

%% Lifts and Drags as function of Angle of Attack and Rotation Speed

% Wing Characteristics (From xflr5)
b     = 0.47 * 2;     % Wingspan
AR    = 7.932;        % Aspect ratio
tr    = 1.37;         % Taper ratio
% S     = 0.111/2;      % Wing Area
% y_MAC = (b / 6) * ((1 + 2 * tr) / (1 + tr)); % Wing location of MAC

% Airfoil Charecteristics
MAC       = 0.119;        % Mean aerodynamic chord
Cl_alpha_0 = 0.0; % From xflr5
Cl_alpha  = 0.1; % From xflr5
Cd0       = 0.01;         % Zero Lift Drag Coefficient
e         = 1.78 * (1 - 0.045 * AR ^.68)-.64; % Oswald's efficiency factor
K         = 1/(pi * AR * e);  % Induced Drag Coefficient
Cl(alpha) = Cl_alpha_0 + Cl_alpha * alpha; % Coefficient of Lift (Symmetric Airfoil)
Cd(alpha) = Cd0 + K * Cl^2;            % Coefficient of Drag

% Wing Coefficients
% CL_alpha_0 = 0.0; % From xflr5
% CL_alpha = 0.0879; % From xflr5
% CL(alpha) = CL_alpha_0 + CL_alpha * alpha;
% LDR(alpha) = 0.001 * alpha^4 - 0.0031 * alpha^3 - ...
%           0.6185 * alpha^2 + 7.1064 * alpha - 5.7227;
% CD(alpha) = CD0 + K * CL^2; 
w_const = [l; S; y_MAC]; % Wing Constants

%% Lift and Drag using formula

v_MAC = (r * y_MAC);    % Velocity of Wing at y_MAC
L1 = (1/2) * CL(alpha1) * S * v_MAC^2;
L2 = (1/2) * CL(alpha2) * S * v_MAC^2;
L3 = (1/2) * CL(alpha3) * S * v_MAC^2;
D1 = (1/2) * CD(alpha1) * S * v_MAC^2;
D2 = (1/2) * CD(alpha2) * S * v_MAC^2;
D3 = (1/2) * CD(alpha3) * S * v_MAC^2;
% D1 = L1 / LDR(alpha1);
% D2 = L2 / LDR(alpha2);
% D3 = L3 / LDR(alpha3);

%% Arm Rotations and Body Rotation

% eta(alpha) = pi/2 - alpha; 
C_B1 = eye(3);              % Arm1 to Body
C_B2 = C(3,deg2rad(120));   % Arm2 to Body
C_B3 = C(3,deg2rad(-120));  % Arm3 to Body
C_1P = C(1,beta(alpha1));   % Prop1 to Arm1
C_2P = C(1,beta(alpha2));   % Prop2 to Arm2
C_3P = C(1,beta(alpha3));   % Prop3 to Arm3
C_BP1 =  C_B1 * C_1P;       % Prop1 to Body
C_BP2 =  C_B2 * C_2P;       % Prop2 to Body
C_BP3 =  C_B3 * C_3P;       % Prop3 to Body

%% Propeller Thrusts

T1_1 = [0; 0; -T1];
T2_2 = [0; 0; -T2];
T3_3 = [0; 0; -T3];
T1_B = C_BP1 * T1_1
T2_B = C_BP2 * T2_2
T3_B = C_BP3 * T3_3
T_c = [T1; T2; T3];

%% Propeller Motor Moments

M(T) = k_M * T; % k_N can be negative if counter-clockwise rotation
M1_1 = M(T1_1);
M2_2 = -M(T2_2);
M3_3 = -M(T3_3);
M1_B = C_BP1 * M1_1
M2_B = C_BP2 * M2_2
M3_B = C_BP3 * M3_3

%% Wing Lifts and Drags (Resulting Moments are Ignored)

L1_B = C_B1 * [0; 0; -L1]
L2_B = C_B2 * [0; 0; -L2]
L3_B = C_B3 * [0; 0; -L3]
D1_B = C_B1 * [0; -D1; 0]
D2_B = C_B2 * [0; -D2; 0]
D3_B = C_B3 * [0; -D3; 0]

%% Wing Blade Motor Moments

W1_1 = [-W1; 0; 0];
W2_2 = [-W2; 0; 0];
W3_3 = [-W3; 0; 0];
W1_B = C_B1 * W1_1
W2_B = C_B2 * W2_2
W3_B = C_B3 * W3_3
W_c = [W1; W2; W3];

%% Arm Length Vectors

l_ = [l; 0; 0];
l1_B = C_B1 * l_
l2_B = C_B2 * l_
l3_B = C_B3 * l_

%% y_MAC Length Vectors

y_ = [y_MAC; 0; 0];
y1_B = C_B1 * y_
y2_B = C_B2 * y_
y3_B = C_B3 * y_

%% Resultant Forces in Body Frame

N_f_g = [0; 0; m * g - R_f];      % Total Gravity Force
B_f_T = T1_B + T2_B + T3_B; % Total Thrust Force
B_f_L = L1_B + L2_B + L3_B; % Total Lift Force
B_f_D = D1_B + D2_B + D3_B; % Total Drag Force
B_f_R = transpose(C_NB) * N_f_g + B_f_T + B_f_L + B_f_D;
B_f_R_symfun = symfun(B_f_R, [T1, T2, T3, alpha1, alpha2, alpha3, ...
                S, y_MAC, phi, theta, r])
B_f_R_hover = simplify(B_f_R_symfun(T1, T2, T3, alpha1, pi/2, pi/2, ...
                S, y_MAC, 0, 0, 0))

%% Resultant Torque in Body Frame

B_tau_T = cross(l1_B, T1_B)...
            + cross(l2_B, T2_B) ...
            + cross(l3_B, T3_B);
B_tau_L = cross(y1_B, L1_B)...
            + cross(y2_B, L2_B) ...
            + cross(y3_B, L3_B);
B_tau_N = M1_B + M2_B + M3_B;
B_tau_M = W1_B + W2_B + W3_B;
B_tau_R = B_tau_T + B_tau_N + B_tau_M;
B_tau_R_symfun = symfun(B_tau_R, [T1, T2, T3, W1, W2, W3, ...
                alpha1, alpha2, alpha3, l])
B_tau_R_hover = simplify(B_tau_R_symfun(T1, T2, T3, 0, 0, 0, alpha1, pi/2, pi/2, l))

%% Equations of Motion

% Primary States  and rate change
rho = [rho_n; rho_e; rho_d]; 
nu = [u; v; w];
Lambda = [phi; theta; psi]; 
omega = [p; q; r]; 
alpha_w = [alpha1; alpha2; alpha3];

% Rate Change States
rho_dot = C_NB * nu
nu_dot = - cross(omega, nu) + (1/m) * B_f_R
Lambda_dot = D * omega
omega_dot =  J \ (B_tau_R - cross(omega, J * omega))
alpha_w_dot = [alpha_dot_1; alpha_dot_2; alpha_dot_3];
alpha_w_ddot = (1 / Iw) * [W1; W2; W3];

% Dynamics State Function
x = [rho; nu; Lambda; omega; alpha_w; alpha_w_dot];
x_dot = simplify([rho_dot; nu_dot; Lambda_dot; omega_dot; alpha_w_dot; alpha_w_ddot]);
u_cB = [B_f_R; B_tau_R];
u_c = [T_c; W_c];

%% Thrust and Moments as function of Motors Inputs

u_c_fun = symfun(u_cB, u_c)
T1 = k1 * delta1;
T2 = k1 * delta2;
T3 = k1 * delta3;
% N1 = k2 * delta1;
% N2 = k2 * delta2;
% N3 = k2 * delta3;
W1 = k3 * sigma1;
W2 = k3 * sigma2;
W3 = k3 * sigma3;
u_c_motor_inputs = u_c_fun(T1, T2, T3, W1, W2, W3);
motor_inputs = [delta1; delta2; delta3; sigma1; sigma2; sigma3];

%% Jacobians (x_dot = A x + B u_c)

A = jacobian(x_dot, x)
B = jacobian(x_dot, u_c)
alloc = jacobian(u_c_motor_inputs, motor_inputs)

%% Saving EOM

sys = [g; m; J_vec; Iw; k_M; R_f];
x_dot_fun = matlabFunction(x_dot, "Vars", {t, x, u_c, sys, w_const})
B_f_R_fun = matlabFunction(B_f_R, "Vars", {t, x, u_c, sys, w_const})
B_tau_R_fun = matlabFunction(B_tau_R, "Vars", {t, x, u_c, sys, w_const})
save("Vortex_EOM2.mat", "x_dot_fun", "B_f_R_fun", "B_tau_R_fun")
