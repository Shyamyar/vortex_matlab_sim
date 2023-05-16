%% Clear
clc;
clear all;
close all;

%% System Parameters
m = 1.59;   % kg
Jn = 232225.1780 * 10^-6;  % kgm^2
Je = 232300.2777 * 10^-6;  % kgm^2
Jd = 454557.3979 * 10^-6;  % kgm^2
g = 9.81;                  % m/sec^s
k_M = 0.1347 / g;       % Slope for Motor Torque vs Thrust (from Thrust Torque measurement experiment)
J_vec= [Jn; Je; Jd];    % Inertia Diagonal Vector
J = [Jn,  0,   0; 
      0, Je,   0; 
      0,  0,  Jd];      % Inertia Tensor
Iw = 645.99 * 10^-6;  % Inertia for Wing at rotation axis
R_f = 0;                % Initial Assumed Reaction force
sys = [g; m; J_vec; Iw; k_M; R_f];
max_T = 0.7 *g;

%% Gain Calculation

% Height gains
zeta_h = 0.707;
a_h3 = 1/m;
h_max = 1;
v_h = sqrt(a_h3*max_T*sqrt(1-zeta_h^2)/h_max);

h_kp = v_h^2/a_h3;
h_kd = 1.1*(2*zeta_h*v_h)/a_h3;
h_ki = 0.3;

% Roll gains
zeta_roll = 0.707;
a_phi3 = 1/Jn;
phi_max = deg2rad(30);
tau_phi_max = 2;
wn_roll = sqrt(a_phi3*tau_phi_max*sqrt(1-zeta_roll^2)/phi_max);

phi_kp = wn_roll^2/a_phi3;
phi_kd = 1.1*(2*zeta_roll*wn_roll)/a_phi3;
phi_ki = 0;

% Pitch gains
zeta_pitch = 0.707;
a_theta3 = 1/Je;
theta_max = deg2rad(30);
tau_theta_max = 2;
wn_pitch = sqrt(a_theta3*tau_theta_max*sqrt(1-zeta_pitch^2)/theta_max);

theta_kp = wn_pitch^2/a_theta3;
theta_kd = 1.1*(2*zeta_pitch*wn_pitch)/a_theta3;
theta_ki = 0;

% Alpha gains
zeta_alpha = 0.707;
a_alpha3 = 1/Iw;
alpha_max = deg2rad(30);
tau_alpha_max = 2;
wn_alpha = sqrt(a_alpha3*tau_alpha_max*sqrt(1-zeta_alpha^2)/alpha_max);

alpha_kp = wn_alpha^2/a_alpha3;
alpha_kd = 1.1*(2*zeta_alpha*wn_alpha)/a_alpha3;
alpha_ki = 2;

%% Lifts and Drags as function of Angle of Attack and Rotation Speed

% Wing Characteristics (From xflr5)
l     = 0.47;         % Arm Length (m)
b     = l * 2;        % Wingspan (m)
AR    = 7.932;        % Aspect ratio
tr    = 1.37;         % Taper ratio
S     = (pi*l^2)/3; %0.116/2;      % Wing Area (m^2)
y_MAC = (b / 6) * ((1 + 2 * tr) / (1 + tr)); % Wing location of MAC

%% Wing Coefficients
w_const = [l; S; y_MAC]; % Wing Constants

%% Initial States
rho0 = [5; 0; -10];
nu0 = [0; 0; 0];
Lambda0 = [0; 0; 0];
omega0 = [0; 0; 0];
alpha0 = deg2rad(90) * ones(3,1); % alpha = 90.55822 for limited_spin, 90 for hover_yaw_drift
alpha_dot_0 = [0; 0; 0];
x0 = [rho0; nu0; Lambda0; omega0; alpha0; alpha_dot_0];
u_c = zeros(6,1);

%% Desired States
a = 5; 
b = 3; 
c = 0; 
T = 10; 
w1 = 2 * pi / T; 
w2 = 0.5 * w1; 
w3 = 1 * w1; 

rhon_d = @(t) a * cos(w2 * t) + 0;     % Desired North Position
rhoe_d = @(t) b * sin(w1 * t) + 0;     % Desired East Position
h_d = @(t) c * cos(w3 * t) + 10;       % Desired Height
w_d = 0;        % Desired vertical speed
alpha_d = 90:45;% AOA fade

%% Load Dependencies
load("Vortex_EOM")

%% Simulation
% Time
t0 = 0;
tend = 50;
t_inc = 0.01;
dt = 0.001;
tspan = 0:t_inc:tend;

% Storage
thist = [];
xhist = [];
u_chist = [];
pids_hist = [];

% Gains
kp_h = h_kp; kd_h = -h_kd; ki_h = h_ki;
kp_roll = phi_kp; kd_roll = phi_kd; ki_roll = phi_ki;
kp_pitch = theta_kp; kd_pitch = theta_kd; ki_pitch = theta_ki;
kp_yaw = 0.073; kd_yaw = 0.057; ki_yaw = 0.0;
kp_alpha = alpha_kp; kd_alpha = alpha_kd; ki_alpha = alpha_ki;

% Modes
% 1. tricopter - regulate thrust, reach constant alpha for no rotation
% 2. vortex - regulate thrust at constant alpha
% 0. null forces and moments - effect of loss of thrust and moments
mode = 1;
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
    alpha = x0(13:15);
    alpha_dot = x0(16:18);
    
    %% Check Ground Reaction, Changing R_f = sys(end)
    sys(end) = 0;
    x0_temp = x0;
    if ground(h)
        x0_temp([3,6]) = 0;
        B_f_R = B_f_R_fun(t0, x0_temp, u_c, sys, w_const);
        N_f_R = C(3,psi) * C(2,theta) * C(1,phi) * B_f_R;
        sys(end) = N_f_R(3);
        if sys(end) >= 0
            x0([3,6]) = 0; % By conservation of momentum (earth has extreme mass)
            h = 0; w = 0;
        end
    end

    %% Changing Modes over time
    if t0<=10
        mode = 1;
    elseif t0<=25
        mode = 1;
    end

    %% Desired States
    T_des = @(alpha) 0; %(m * g) / 3 - (S*CL(alpha)*r^2*y_MAC^2)/2; % Desired Thrust per prop for hover
    if mode == 1
        phi_d = 0.075 * (rhoe_d(tspan(i)) - rhoe) - 0.15 * v; p_d = 0; % roll control
        theta_d = - 0.075 * (rhon_d(tspan(i)) - rhon) + 0.15 * u; q_d = 0; % pitch control
%         phi_d = 0; p_d = 0; % roll control
%         theta_d = 0; q_d = 0; % pitch control
        psi_d = 0; r_d = 0; % yaw control
    elseif mode == 2
        phi_d = 0.075 * (rhoe_d(tspan(i)) - rhoe) - 0.15 * v; p_d = 0; % roll control
        theta_d = - 0.075 * (rhon_d(tspan(i)) - rhon) + 0.15 * u; q_d = 0; % pitch control
%         phi_d = 0; p_d = 0; % roll control
%         theta_d = 0; q_d = 0; % pitch control
        psi_d = psi; r_d = r; % no yaw control
    end
    
    %% State Errros
    e_h = h_d(tspan(i)) - h;
    e_roll = phi_d - phi;
    e_pitch = theta_d - theta;
    e_yaw = psi_d - psi;
    udot_prev = u;
    vdot_prev = v;
    if i == 1
        e_h_int = 0;
        e_roll_int = 0;
        e_pitch_int = 0;
        e_yaw_int = 0;
        e_alpha_int = 0;
    end

    %% PID Control with desired path
    if mode == 1
        k_h_mode = 1;
        k_att_mode = 1;
    elseif mode == 2
        k_h_mode = 50;
        k_att_mode = 0.01;
    end
    h_pid = k_h_mode * ( kp_h * e_h + kd_h * (w_d - w) + ki_h * e_h_int);
    roll_pid = k_att_mode * (kp_roll * e_roll + kd_roll * (p_d - p) + ki_roll * e_roll_int);
    pitch_pid = k_att_mode * (kp_pitch * e_pitch + kd_pitch * (q_d - q) + ki_pitch * e_pitch_int);
    yaw_pid1 = kp_yaw * e_yaw + kd_yaw * (r_d - r) + ki_yaw * e_yaw_int;
    yaw_pid2 = 0.5 * e_yaw + 0.2 * (r_d - r) + 10 * e_yaw_int;
    
    %% Alpha for yaw control
    if mode == 1
        alpha_c = alpha0 - yaw_pid1 * ones(3,1) + 0.0097;
        e_alpha = alpha_c - alpha;
        T1 = min(max(T_des(alpha(1)) + h_pid - pitch_pid, 0), max_T);
        T2 = min(max(T_des(alpha(2)) + h_pid + 0.5 * pitch_pid + roll_pid, 0), max_T);
        T3 = min(max(T_des(alpha(3)) + h_pid + 0.5 * pitch_pid - roll_pid, 0), max_T);
        W = kp_alpha * (e_alpha) + kd_alpha * (0 - alpha_dot) + ki_alpha * e_alpha_int;

    elseif mode == 2
        alpha_c = max(alpha - 0.05, deg2rad(45) * ones(3,1));
        e_alpha = alpha_c - alpha;
        psi_lag = psi - 0.0;
        T1 = min(max(T_des(alpha(1)) + h_pid - pitch_pid * cos(psi_lag), 0), max_T);
        T2 = min(max(T_des(alpha(2)) + h_pid + 0.5 * pitch_pid * cos(psi_lag) + ...
                roll_pid * cos(psi_lag), 0), max_T);
        T3 = min(max(T_des(alpha(3)) + h_pid + 0.5 * pitch_pid * cos(psi) - ...
                roll_pid * cos(psi_lag), 0), max_T);
        W = kp_alpha * (e_alpha) + kd_alpha * (0 - alpha_dot) + ki_alpha * e_alpha_int;
    else
        T1 = 0; T2 = 0; T3 = 0;
        W = zeros(3,1);
    
    end
    
    %% Integral part calculation
    e_h_int = e_h_int + e_h * t_inc;
    e_roll_int = e_roll_int + e_roll * t_inc;
    e_pitch_int = e_pitch_int + e_pitch * t_inc;
    e_yaw_int = e_yaw_int + e_yaw * t_inc;
    e_alpha_int = e_alpha_int + e_alpha * t_inc;

    %% Consolidate Inputs
    T_c = [T1; T2; T3];
    W_c = W;
    u_c = [T_c; W_c];
    pids = [h_pid; pitch_pid; roll_pid];

    %% Dynamics/Plant
    [t, x] = ode45(x_dot_fun, t0:dt:t0+t_inc, x0', [], u_c, sys, w_const);

    %% Storage
    x(:, 7:9) = wrapToPi(x(:, 7:9)); % wrapping Lambda between -pi to pi
    x(:, 13:15) = wrapTo2Pi(x(:, 13:15)); % wrapping alpha between 0 to 2pi
    thist = [thist; t(1:end-1)];
    xhist = [xhist; x(1:end-1, :)];
    u_chist = [u_chist, u_c];
    pids_hist = [pids_hist, pids];
    x0 = x(end, :)';
    t0 = t0 + t_inc;
    
    if ground(h) && t0 > 0.5
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
    hold on
    if i == 1
        plot(tspan, rhon_d(tspan))
    elseif i == 2
        plot(tspan, rhoe_d(tspan))
    elseif i == 3
        plot(tspan, h_d(tspan))
    end
    xlabel("time (sec)");
    ylabel(x_name(i),'Interpreter','latex')
    grid on
end
return

%% Control Plots
fig2 = figure(2);
T_leg = ["T1", "T2", "T3"];
plot(tspan,u_chist(1:3,:)'/g)
% yyaxis("right")
% hold on
% plot(thist, xhist(:,[9,13]))
% yyaxis("left")
legend(T_leg)

fig3 = figure(3);
pids_leg = ["h_pid", "pitch_pid", "roll_pid"];
plot(tspan,pids_hist(2:3,:)')
legend(pids_leg(2:3))


%% Animate
VIDEO = false;
video_title = "..\Videos\hover_drift";
fig3 = animate_system_seq(xhist, thist, VIDEO, video_title, 0.01);

%% Functions
function ground = ground(h)
ground = true;    
    if h>0
        ground = false;
    end
end