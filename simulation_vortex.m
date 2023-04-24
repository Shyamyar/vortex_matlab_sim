%% Clear
clc;
clear all;
close all;


%% System Parameters
m = 1.59;   % kg
Jx = 525154.5418 * 10^-6;  % kgm^2
Jy = 640185.2986 * 10^-6;  % kgm^2
Jz = 127840.6707 * 10^-6;  % kgm^2
g = 9.81;           % m/sec^s
k_M = 0.1347 / g;     % Slope for Motor Torque vs Thrust (from Thrust Torque measurement experiment)
J_vec= [Jx; Jy; Jz];% Inertia Diagonal Vector
J = [Jx,  0,   0; 
      0, Jy,   0; 
      0,  0,  Jz];  % Inertia Tensor
Iw = 673.3897 * 10^-6;          % Inertia for Wing at rotation axis
sys = [g; m; J_vec; Iw; k_M];

%% Lifts and Drags as function of Angle of Attack and Rotation Speed

% Wing Characteristics (From xflr5)
l     = 0.47;         % Arm Length (m)
b     = l * 2;        % Wingspan (m)
AR    = 7.932;        % Aspect ratio
tr    = 1.37;         % Taper ratio
S     = 0.116/2;      % Wing Area (m^2)
y_MAC = (b / 6) * ((1 + 2 * tr) / (1 + tr)); % Wing location of MAC

%% Wing Coefficients
% CL = CL_fun(alpha);
% CD = CD_fun(alpha);
% LDR = LDR_fun(alpha);
w_const = [l; S; y_MAC]; % Wing Constants

%% Initial States
rho0 = [0; 0; 0];
nu0 = [0; 0; 0];
Lambda0 = [0; 0; 0];
omega0 = [0; 0; 0];
alpha0 = deg2rad(90) * ones(3,1); % alpha1 = 89.8 for limited_spin, 90 for hover_yaw_drift
alpha_dot_0 = [0; 0; 0];
x0 = [rho0; nu0; Lambda0; omega0; alpha0; alpha_dot_0];

%% Desired States
rhon_d = 0;     % Desired North Position
rhoe_d = 0;     % Desired East Position
h_d = 10;       % Desired Height
w_d = 0;        % Desired vertical speed

%% Load Dependencies
load("Vortex_EOM")

%% Simulation
% Time
t0 = 0;
tend = 10;
t_inc = 0.01;
dt = 0.001;
tspan = 0:t_inc:tend;

% Storage
thist = [];
xhist = [];
Thist = [];

% Gains
kp_h = 3.5; kd_h = -2; ki_h = 100;
kp_roll = 0.2; kd_roll = 0.05; ki_roll = 0.345;
kp_pitch = 0.2; kd_pitch = 0.05; ki_pitch = 0.345;
kp_yaw = 0.3; kd_yaw = 0.5; ki_yaw = 0.345;

% Modes
tricopter = true;
vortex = ~tricopter;
PID = true;

for i = 1:length(tspan)
    %% Named States
    rhon = x0(1);
    rhoe = x0(2);
    h = -x0(3);
    u = x0(4);
    v = x0(5);
    w = x0(6);    
    phi = x0(7);
    theta = x0(8);
    psi = x0(9);    
    p = x0(10);
    q = x0(11);
    r = x0(12);
    alpha1 = x0(13);
    alpha2 = x0(14);
    alpha3 = x0(15);
    alpha_dot1 = x0(16);
    alpha_dot2 = x0(17);
    alpha_dot3 = x0(18);

    %% Control Inputs
    T_des = @(alpha) 0; %(m * g) / 3 - (S*CL(alpha)*r^2*y_MAC^2)/2; % Desired Thrust per prop for hover
    if tricopter
        phi_d = 0; p_d = 0; % roll control
        theta_d = 0; q_d = 0; % pitch control
        psi_d = 0; r_d = 0; % yaw control
    end
    if vortex
        phi_d = 0; p_d = 0; % roll control
        theta_d = 0; q_d = 0; % pitch control
        psi_d = psi; r_d = r; % no yaw control
    end
    
    %% State Errros
    e_h = h - h_d;
    e_roll = phi - phi_d;
    e_pitch = theta - theta_d;
    e_yaw = psi - psi_d;
    if i == 1
        e_h_prev = 0;
        e_roll_prev = 0;
        e_pitch_prev = 0;
        e_yaw_prev = 0;
    end
    e_h_int = e_h_prev + e_h * t_inc;
    e_roll_int = e_roll_prev + e_roll * t_inc;
    e_pitch_int = e_pitch_prev + e_pitch * t_inc;
    e_yaw_int = e_yaw_prev + e_yaw * t_inc;

    %% PID Control with desired path
    h_pid = 12.5 * (- kp_h * e_h - kd_h * (w - w_d) - ki_h * e_h_int);
    roll_pid = 0.1 * (- kp_roll * e_roll - kd_roll * (p - p_d) - ki_roll * e_roll_int);
    pitch_pid = 0.1 * (- kp_pitch * e_pitch - kd_pitch * (q - q_d) - ki_pitch * e_pitch_int);
    yaw_pid1 = 0.1 * (- kp_yaw * (psi - psi_d) - kd_yaw * (r - r_d) - ki_yaw * e_yaw_int);
    yaw_pid2 = - 0.05 * (psi - psi_d) - 0.09 * (r - r_d) - 10 * e_yaw_int;
    
    %% Alpha for yaw control
    if tricopter
%         yaw_pid = 0; % no yaw control
        alpha_c = alpha0 - yaw_pid1;
        T1 = min(T_des(alpha1) + h_pid - pitch_pid, 0.7*g);
        T2 = min(T_des(alpha2) + h_pid + pitch_pid - roll_pid, 0.7*g);
        T3 = min(T_des(alpha3) + h_pid + pitch_pid + roll_pid, 0.7*g);
%         alpha_c = deg2rad(90.55822) * ones(3,1);
    end
    if vortex
        alpha_c = deg2rad(15) * ones(3,1);
        T1 = min(T_des(alpha1) + h_pid - pitch_pid * cos(psi), 0.7*g);
        T2 = min(T_des(alpha2) + h_pid + pitch_pid * cos(psi) - ...
                roll_pid * cos(psi), 0.7*g);
        T3 = min(T_des(alpha3) + h_pid + pitch_pid * cos(psi) + ...
                roll_pid * cos(psi), 0.7*g);
    end

    if PID
        kp_W = 0.05; % 0.03
        kd_W = 0.010; % 0.0
        W1 = - kp_W * (alpha1 - alpha_c(1)) - kd_W * (alpha_dot1 - 0);
        W2 = - kp_W * (alpha2 - alpha_c(2)) - kd_W * (alpha_dot2 - 0);
        W3 = - kp_W * (alpha3 - alpha_c(3)) - kd_W * (alpha_dot3 - 0);
%         W1 = -yaw_pid2;
%         W2 = -yaw_pid2;
%         W3 = -yaw_pid2;
        
    else
        %% Direct Inputs
        T1 = T_des(alpha1); % Tper = 2 for limited_spin, 1 for hover_yaw_drift
        T2 = T_des(alpha2) - roll_pid * cos(psi);
        T3 = T_des(alpha3) + roll_pid * cos(psi);
        W1 = 0;
        W2 = 0;
        W3 = 0;
    end
    
    %% Store previous state errors
    e_h_prev = e_h * t_inc;
    e_roll_prev = e_roll * t_inc;
    e_pitch_prev = e_pitch * t_inc;
    e_yaw_prev = e_yaw * t_inc;

    %% Consolidate Inputs
    T_c = [T1; T2; T3];
    W_c = [W1; W2; W3];
    u_c = [T_c; W_c];

    %% Dynamics/Plant
    [t, x] = ode45(x_dot_fun, t0:dt:t0+t_inc, x0', [], u_c, sys, w_const);

    %% Storage
    x(:, 7:9) = wrapToPi(x(:, 7:9)); % wrapping Lambda between -pi to pi
    x(:, 13:15) = wrapTo2Pi(x(:, 13:15)); % wrapping alpha between 0 to 2pi
    thist = [thist; t(1:end-1)];
    xhist = [xhist; x(1:end-1, :)];
    Thist = [Thist, repmat([T1; T2; T3], 1, size(t,1)-1)];
    x0 = x(end, :)';
    t0 = t0 + t_inc;
    
    if h > 20 || h < -20
        disp(tspan(i));
        break;
    end
end

%% State Plots
x_name = [...
            "$\rho_n (m)$"      "$\rho_e (m)$"      "h (m)" ...
            "u (m/s)"           "v (m/s)"           "w (m/s)" ...
            "$\phi (rad)$"      "$\theta (rad)$"    "$\psi (rad)$" ...
            "p (rad/s)"         "q (rad/s)"         "r (rad/s)" ...
            "$\alpha_1 (rad)$"  "$\alpha_2 (rad)$"  "$\alpha_3 (rad)$" ...
            "$\dot\alpha_1  (rad/s)$"    "$\dot\alpha_2 (rad/s)$"    "$\dot\alpha_3 (rad/s)$"];

fig1 = figure("WindowState","maximized");
xplot = xhist;
xplot(:, 3) = - xhist(:, 3); % wrapping Lambda between -pi to pi
xplot(:, 7:9) = wrapToPi(xhist(:, 7:9)); % wrapping Lambda between -pi to pi
xplot(:, 13:15) = wrapTo2Pi(xhist(:, 13:15)); % wrapping alpha between 0 to 2pi

for i = 1:18
    subplot(6,3,i)
    plot(thist, xplot(:, i))
    xlabel("time (sec)");
    ylabel(x_name(i),'Interpreter','latex')
end
return
%% Control Plots
fig2 = figure(2);
T_leg = ["T1", "T2", "T3"];
plot(thist,Thist'/g)
% yyaxis("right")
% hold on
% plot(thist, xhist(:,[9,13]))
% yyaxis("left")
legend(T_leg)

%% Animate
VIDEO = false;
video_title = "hover_drift";
fig3 = animate_system_seq(xhist, thist, VIDEO, video_title, 0.01);