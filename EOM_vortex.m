%% State Variables
syms t m Jn Je Jd g Iw
syms rho_n rho_e rho_d u v w phi theta psi p q r 
syms kp_phi kp_theta kp_psi kd_phi kd_theta kd_psi

%% Symbols for Forces, Moments, Control Inputs
syms alpha % AOA
syms alpha1 alpha2 alpha3  % Wing position w.r.t vertical
syms alpha_dot_1 alpha_dot_2 alpha_dot_3 % rate of change of alpha
% syms Cl_alpha_0 Cl_alpha Cd0 K % Lift and Drag Coefficient elements
syms CL_alpha_0  CL_alpha % Zero lift Intercept and Slope 
syms l S y_MAC MAC % arm length to tip motor, Wing area, Wing location of MAC 
% syms Cl1 Cl2 Cl3 Cd1 Cd2 Cd3 % Individual coefficients of Lift and Drag
syms LDR(alpha) % Lift to Drag ratio
syms Cl(alpha) Cd(alpha) % Airfoil coefficient of lift and drag
syms CL(alpha) CD(alpha) % Wing coefficient of lift and drag
syms CLx(alphaeff)
syms CDx(alphaeff)
syms Cmx(alphaeff)
syms alpha_eff(vec1,vec2,vec3)
syms T1 T2 T3 T % Thrust for each propeller
syms M(T) k_M % Moments for each propeller rotation proportional to Thrust
syms W1 W2 W3 W % Moments for each wing rotation
syms delta1 delta2 delta3 % tip motor inputs
syms sigma1 sigma2 sigma3 % servo motor inputs
syms k1 k2 k3 % Gains mapping to Thrust, Motor Moment, and Wing Moment
syms R_f % Reaction force in vertical direction

%% System Parameters

J_vec= [Jn; Je; Jd]; 
J = [Jn,  0,   0; 
      0, Je,   0; 
      0,  0,  Jd]; 
C_NB = C(3, psi) * C(2, theta) * C(1, phi); %DCM to convert from body to inertial
D = D_matrix(phi, theta);

%% Lifts and Drags as function of Angle of Attack and Rotation Speed

% % Wing Characteristics (From xflr5)
% b     = 0.47 * 2;     % Wingspan
% AR    = 7.932;        % Aspect ratio
% tr    = 1.37;         % Taper ratio
% S     = 0.111/2;      % Wing Area
% y_MAC = (b / 6) * ((1 + 2 * tr) / (1 + tr)); % Wing location of MAC
% MAC       = 0.124;        % Mean aerodynamic chord
% 
% % Airfoil Charecteristics
% Cl_alpha_0 = 0.0; % From xflr5
% Cl_alpha  = 0.1; % From xflr5
% Cd0       = 0.01;         % Zero Lift Drag Coefficient
% e         = 1.78 * (1 - 0.045 * AR ^.68)-.64; % Oswald's efficiency factor
% K         = 1/(pi * AR * e);  % Induced Drag Coefficient
% Cl(alpha) = Cl_alpha_0 + Cl_alpha * alpha; % Coefficient of Lift (Symmetric Airfoil)
% Cd(alpha) = Cd0 + K * Cl^2;            % Coefficient of Drag
% 
% % Wing Coefficients
% CL_alpha_0 = 0.0; % From xflr5
% CL_alpha = 0.0879; % From xflr5
% CL(alpha) = CL_alpha_0 + CL_alpha * alpha;
% LDR(alpha) = 0.001 * alpha^4 - 0.0031 * alpha^3 - ...
%           0.6185 * alpha^2 + 7.1064 * alpha - 5.7227;
% CD(alpha) = CD0 + K * CL^2; 

%% Arm Rotations and Body Rotation

C_flip = [1 0 0;
          0 1 0;
          0 0 -1];
C_B1 = C(3,pi);      % Arm1 to Body
C_B1([2,3,4]) = 0;
C_B2 = C(3,-pi/3);   % Arm2 to Body
C_B3 = C(3,pi/3);   % Arm3 to Body
C_1P = C(1,alpha1) * C_flip;   % Prop1 to Arm1
C_2P = C(1,alpha2) * C_flip;   % Prop2 to Arm2
C_3P = C(1,alpha3) * C_flip;   % Prop3 to Arm3
C_B1P =  C_B1 * C_1P;       % Prop1 to Body
C_B2P =  C_B2 * C_2P;       % Prop2 to Body
C_B3P =  C_B3 * C_3P;       % Prop3 to Body

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

%% MAC Length Vectors

mac_ = [0; 0; -MAC/2]; % Assuming pivot axis is at half the MAC
mac_1P = C_1P * mac_;
mac_2P = C_2P * mac_;
mac_3P = C_3P * mac_;

%% Primary States
rho = [rho_n; rho_e; rho_d]; 
nu = [u; v; w];
Lambda = [phi; theta; psi]; 
omega = [p; q; r]; 
alpha_w = [alpha1; alpha2; alpha3];
alpha_w_dot = [alpha_dot_1; alpha_dot_2; alpha_dot_3];
x = [rho; nu; Lambda; omega; alpha_w; alpha_w_dot];

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
T_c = [T1; T2; T3]

%% Propeller Motor Moments

M(T) = k_M * T; % k_M can be negative if counter-clockwise rotation
M1_1 = M(T1_1);
M2_2 = -M(T2_2);
M3_3 = -M(T3_3);
M1_B = C_B1P * M1_1
M2_B = C_B2P * M2_2
M3_B = C_B3P * M3_3

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

% B_f_R_symfun = symfun(B_f_R, [T1, T2, T3, alpha1, alpha2, alpha3, ...
%                 S, y_MAC, phi, theta, r])
% B_f_R_hover = simplify(B_f_R_symfun(T, T, T, pi/2, pi/2, pi/2, ...
%                 S, y_MAC, 0, 0, 0))
% B_f_R_vortex = simplify(B_f_R_symfun(T, T, T, pi/4, pi/4, pi/4, ...
%                 S, y_MAC, phi, theta, r))

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

% B_tau_R_symfun = symfun(B_tau_R, [T1, T2, T3, W1, W2, W3, ...
%                 alpha1, alpha2, alpha3, l])
% B_tau_R_hover = simplify(B_tau_R_symfun(T, T, T, 0, 0, 0, alpha1, pi/2, pi/2, l))

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
u_cB = [B_f_R; B_tau_R];
u_c = [T_c; W_c];

%% Thrust and Moments as function of Motors Inputs

% u_c_fun = symfun(u_cB, u_c)
% T1 = k1 * delta1;
% T2 = k1 * delta2;
% T3 = k1 * delta3;
% N1 = k2 * delta1;
% N2 = k2 * delta2;
% N3 = k2 * delta3;
% W1 = k3 * sigma1;
% W2 = k3 * sigma2;
% W3 = k3 * sigma3;
% u_c_motor_inputs = u_c_fun(T1, T2, T3, W1, W2, W3);
% motor_inputs = [delta1; delta2; delta3; sigma1; sigma2; sigma3];

%% Jacobians (x_dot = A x + B u_c)

% A = jacobian(x_dot, x)
% B = jacobian(x_dot, u_c)
% alloc = jacobian(u_c_motor_inputs, motor_inputs)

%% Saving EOM

sys = [g; m; J_vec; Iw; k_M; R_f];
w_const = [l; S; y_MAC; MAC]; % Wing Constants
x_dot_fun = matlabFunction(x_dot, "Vars", {t, x, u_c, sys, w_const})
B_f_R_fun = matlabFunction(B_f_R, "Vars", {t, x, u_c, sys, w_const})
B_tau_R_fun = matlabFunction(B_tau_R, "Vars", {t, x, u_c, sys, w_const})
save("Vortex_EOM.mat", "x_dot_fun", "B_f_R_fun", "B_tau_R_fun")
