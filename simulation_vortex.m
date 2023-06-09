%% Clear
clc;
clear;
close all;

%% Drone Body Parameters
% Wingold_properties;
% Wing550_properties;
Wing350_properties;

%% System Parameters
g = 9.81;                  % m/sec^s
k_M = 0.1347 / g;       % Slope for Motor Torque vs Thrust (from Thrust Torque measurement experiment)
J_vec= [Jn; Je; Jd];    % Inertia Diagonal Vector
J = [Jn,  0,   0; 
      0, Je,   0; 
      0,  0,  Jd];      % Inertia Tensor
R_f = m*g;              % Initial Assumed Reaction force
max_T = 0.7 * g;
max_W = 0.6278;         % Hitec HS-7245MH
max_vel = 2;
tot_mAh = 6000;         % Battery Capacity
tot_V = 22.2;           % Battery Operating Voltage
tot_Wh = tot_V * tot_mAh * 1e-3;

%% Wing Coefficients
w_const = [l; S; y_MAC; MAC]; % Wing Constants

%% Load Dependencies
load("Vortex_EOM")

%% Storage Variables w.r.t AOA
r_alpha_hist = [];
T_alpha_hist = [];
W_alpha_hist = [];
u_chist_mode = [];
tot_flight_time_mode = [];
cum_Wh_mode = [];

%% Initial States
rho0 = [0; 0; 0];
nu0 = [0; 0; 0];
Lambda0 = [0; 0; 0];
omega0 = [0; 0; 0];
alpha0 = deg2rad([0;0;0]); % alpha = 90.55822 for limited_spin, 90 for hover_yaw_drift
alpha_dot_0 = [0; 0; 0];

% rho0 = [0; 0; -10];
% Lambda0 = [0; -pi/2.5; 0];
% C_NB = C(3,Lambda0(3)) * C(2,Lambda0(2)) * C(1,Lambda0(1));
% nu0 = transpose(C_NB) * [10;0;0]; % [5;0;0] in inertial frame
% omega0 = [0; 0; 0];
% alpha0 = deg2rad([0;-5;5]); % alpha = 90.55822 for limited_spin, 90 for hover_yaw_drift
% alpha_dot_0 = [0; 0; 0];

%% Desired States
ns = 0; 
es = 0; 
hs = 0; 
T = 10; 
w1 = 2 * pi / T; 
w2 = 0.5 * w1; 
w3 = 1 * w1; 
limit = @(val,lim) min(max(val,-lim),lim);

rhon_d = @(t,n0,ns) ns * cos(w2 * t) + n0;     % Desired North Position
rhoe_d = @(t,e0,es) es * sin(w1 * t) + e0;     % Desired East Position
h_d = @(t,h0,hs) hs * cos(w3 * t) + h0;       % Desired Height

%% Simulation
alphas = 45;
modes = 1;
t_inc = 0.01;
tend = 50;
tspan = 0:t_inc:tend;

for mode_id = 1:numel(modes)

for alpha_des_id = 1:numel(alphas)
tic;
alpha_des = alphas(alpha_des_id);

%% Intial States, Controls and System
x0 = [rho0; nu0; Lambda0; omega0; alpha0; alpha_dot_0];
u_c = [4*ones(3,1);zeros(3,1)];
n0 = 0;
e0 = 0;
h0 = 10;
way_point = [rhon_d(0,n0,ns); rhoe_d(0,e0,es); h_d(0,h0,hs)];
sys = [g; m; J_vec; Iw; k_M; R_f];

%% Time Parameters
t0 = 0;
dt = 0.001;

%% State Storage Initialization
thist = [];
xhist = [];
u_chist = [];
pids_hist = [];
deshist = [];

%% Integral errors Initialization
e_n_int = 0;
e_e_int = 0;
e_h_int = 0;
e_w_int = 0;
e_roll_int = 0;
e_pitch_int = 0;
e_yaw_int = 0;
e_alpha_int = 0;

%% Other Initializations
psi_l_prev = 0;

%% Start Mode
% Modes
% 1. tri-copter - regulate thrust, reach constant alpha for no rotation
% 2. windmill - regulate thrust at constant alpha
% 3. fixed-wing - rotate at certain pitch angle and increase alpha for lift
% 0. null or constant forces and moments - effect of loss of thrust and moments

start_mode = modes(mode_id); % Start mode
mode = start_mode;
modefixed = false;

for i = 1:length(tspan)
    disp(tspan(i))
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
    x0_temp = x0;
    sys(end) = 0;
    if ground(h)
        x0_temp([3,6]) = 0;
        B_f_R = B_f_R_fun(t0, x0_temp, u_c, sys, w_const);
        N_f_R = C(3,psi) * C(2,theta) * C(1,phi) * B_f_R;
        sys(end) = N_f_R(3);
        if sys(end) >= 0
            x0([3,6]) = 0; % By conservation of momentum (earth has extreme mass)
            h = 0;
        end
    end

    %% Changing Modes over time
    mode_prev = mode;
    if tspan(i) >= 20
        mode = 2;
        t_trans = 2; % sec
        r_trans = r;
    end
    
    if mode ~= mode_prev || i == 1 % Get gains everytime mode changes
        [kp, kd, ki] = run_gains(m, Jn, Je, Jd, Iw, max_T, max_W, max_vel, mode, alpha_des);
    end

    %% Changing Waypoints
    if tspan(i) >= 10
        % ns = 4; 
        % es = 3; 
        % hs = 0.5;
        % n0 = 3;
        % e0 = 0;
        % h0 = 10;
        % way_point = [rhon_d(tspan(i),n0,ns); rhoe_d(tspan(i),e0,es); h_d(tspan(i),h0,hs)];
    end
    
    %% Lookahead Point
    ldist = norm([1.4,1.4]);
    trans_t = way_point - [rhon; rhoe; h];       % relative vector to target
    d_slant_t = norm(trans_t);           % distance to target
    unit_trans_t = trans_t ./ d_slant_t;            % unit vector signifying direction towards target
    if any(isnan(unit_trans_t))
        unit_trans_t = zeros(3,1);
    end
    if norm(trans_t(1:2)) > 0.1 && mode == 1
        psi_l = atan2(unit_trans_t(2), unit_trans_t(1));    % desired course angle towards target
        if psi_l_prev - psi_l < -pi
            psi_l = psi_l - 2*pi;
        elseif psi_l_prev - psi_l > pi
            psi_l = psi_l + 2*pi;
        end
    else
        psi_l = psi;
    end
    psi_l_prev = psi_l;
    lookahead_dist = min(d_slant_t, ldist);          % Minimum of threshold lookahead and nearest target point
    del_pos_l = lookahead_dist .* unit_trans_t;     % Transition Target Position w.r.t current position
    l_des = [rhon; rhoe; h] + del_pos_l;     % Transition Target Position in inertial frame
    deshist = [deshist, [way_point; psi_l]];

    %% Desired States
    if mode == 1
        des_pos_l_ = C(3,psi).' * del_pos_l;
        e_n = des_pos_l_(1);
        e_e = des_pos_l_(2);
        e_h = des_pos_l_(3);
        u_d = limit(0.2 * e_n, max_vel);
        v_d = limit(0.2 * e_e, max_vel);
        w_d = -limit(1 * e_h, max_vel);
        e_u = u_d - u;
        e_v = v_d - v;
        e_w = w_d - w;
        psi_d = psi_l;
        r_d = 0; % yaw control
        phi_d = limit(kp.e * e_e + kd.e * e_v + ki.e * e_e_int, deg2rad(10)); 
        p_d = 0; % roll control
        theta_d = limit(- kp.n * e_n - kd.n * e_u  - ki.n * e_n_int, deg2rad(10)); 
        q_d = 0; % pitch control
    elseif mode == 2 || mode == 0
        des_pos_l_ = C(3,0).' * del_pos_l;
        e_n = des_pos_l_(1);
        e_e = des_pos_l_(2);
        e_h = des_pos_l_(3);
        u_d = limit(0.2 * e_n, max_vel);
        v_d = limit(0.2 * e_e, max_vel);
        w_d = -limit(1 * e_h, max_vel);
        e_u = u_d - u;
        e_v = v_d - v;
        e_w = w_d - w;
        psi_d = psi_l; r_d = r; % no yaw control
        phi_d = limit(kp.e * e_e + kd.e * e_v + ki.e * e_e_int, deg2rad(30)); 
        p_d = 0; % roll control
        theta_d = limit(- kp.n * e_n - kd.n * e_u  - ki.n * e_n_int, deg2rad(30));
        q_d = 0; % pitch control
    elseif mode == 3
        des_pos_l_ = C(3,0).' * del_pos_l;
        e_n = 0;
        e_e = 0;
        e_h = des_pos_l_(3);
        theta_d = max(theta - 0.2, deg2rad(-65)); 
        q_d = 0; % pitch control
        vel_d = 10;
        u_d = vel_d * cos(theta_d);
        v_d = 0;
        w_d = vel_d * sin(theta_d);
        e_u = u_d - u;
        e_v = v_d - v;
        e_w = w_d - w;
        phi_d = limit(kp.e * e_e + kd.e * e_v + ki.e * e_e_int, deg2rad(30)); 
        p_d = 0; % roll control
        psi_d = 0; r_d = 0; % yaw control
    end
    
    %% Attitude State Errors
    e_roll = phi_d - phi;
    e_pitch = theta_d - theta;
    e_yaw = psi_d - psi;
    uprev = u;
    vprev = v;
    wprev = w;

    %% PID Control with desired path
    h_pid = kp.h * e_h + kd.h * e_w + ki.h * e_h_int;
    roll_pid = kp.roll * e_roll + kd.roll * (p_d - p) + ki.roll * e_roll_int;
    pitch_pid = kp.pitch * e_pitch + kd.pitch * (q_d - q) + ki.pitch * e_pitch_int;
    yaw_pid = kp.yaw * e_yaw + kd.yaw * (r_d - r) + ki.yaw * e_yaw_int;
    
    %% Alpha for yaw control
    if mode == 1
        T1 = min(max(h_pid - pitch_pid, 0), max_T);
        T2 = min(max(h_pid + 0.5 * pitch_pid + roll_pid, 0), max_T);
        T3 = min(max(h_pid + 0.5 * pitch_pid - roll_pid, 0), max_T);

        alpha_c =  yaw_pid * ones(3,1) + 0.0097;
        % if (mode_prev == 2 || t_trans > 0) && r > 0
        %     t_trans = t_trans - 0.1;
        %     alpha_c = -(deg2rad(2)/r_trans) * r;
        % end
        e_alpha = alpha_c - alpha;
        W_c = limit(kp.alpha * (e_alpha) + kd.alpha * (0 - alpha_dot) ...
            + ki.alpha * e_alpha_int, max_W);

    elseif mode == 2
        psi_lag = psi + deg2rad(90); 
        % if tspan(i)>=15
        %     roll_pid = 0;
        %     pitch_pid = -7;
        % end
        T1 = min(max(h_pid - pitch_pid * cos(psi_lag) + ...
                roll_pid * sin(psi_lag) , 0), max_T);
        T2 = min(max(h_pid - pitch_pid * cos(psi_lag + 2 * pi/3) + ...
                roll_pid * sin(psi_lag + 2 * pi/3) , 0), max_T);
        T3 = min(max(h_pid - pitch_pid * cos(psi_lag - 2 * pi/3) + ...
                roll_pid * sin(psi_lag - 2 * pi/3), 0), max_T);

        alpha_c = deg2rad(alpha_des) * ones(3,1);
        e_alpha = alpha_c - alpha;
        W_c = limit(kp.alpha * (e_alpha) + kd.alpha * (0 - alpha_dot) ...
            + ki.alpha * e_alpha_int, max_W);
    
    elseif mode == 3
        h_pid = 0.9 * kp.h * e_h + 0.77 * kd.h * e_w;
        T1 = min(max(h_pid - pitch_pid, 0), max_T);
        T2 = min(max(h_pid + 0.5 * pitch_pid + roll_pid, 0), max_T);
        T3 = min(max(h_pid + 0.5 * pitch_pid - roll_pid, 0), max_T);

        alpha_c =  yaw_pid * ones(3,1) + 0.0097;
        aoa = deg2rad(5);
        alpha_c = alpha_c + [0; min(max(alpha(2) - aoa/100, -aoa),0);...
            max(min(alpha(3) + aoa/100, aoa),0)];
        e_alpha = alpha_c - alpha;
        W_c = limit(kp.alpha * (e_alpha) + kd.alpha * (0 - alpha_dot) ...
            + ki.alpha * e_alpha_int, max_W);

    elseif mode == 0
        T1 = 2; T2 = 2; T3 = 2;
        alpha_c = min(alpha + 0.05, deg2rad(alpha_des) * ones(3,1)); % AOA fade
        e_alpha = alpha_c - alpha;
        W_c = limit(kp.alpha * (e_alpha) + kd.alpha * (0 - alpha_dot) ...
            + ki.alpha * e_alpha_int, max_W);
    
    end

    %% Differential part calculation
    udot = (u - uprev)/t_inc;
    vdot = (v - vprev)/t_inc;
    wdot = (w - wprev)/t_inc;

    %% Integral part calculation
    e_n_int = e_n_int + e_n * t_inc;
    e_e_int = e_e_int + e_e * t_inc;
    e_h_int = e_h_int + e_h * t_inc;
    e_w_int = e_w_int + e_w * t_inc;
    e_roll_int = e_roll_int + e_roll * t_inc;
    e_pitch_int = e_pitch_int + e_pitch * t_inc;
    e_yaw_int = e_yaw_int + e_yaw * t_inc;
    e_alpha_int = e_alpha_int + e_alpha * t_inc;

    %% Consolidate Inputs
    T_c = [T1; T2; T3];
    u_c = [T_c; W_c];
    pids = [h_pid; pitch_pid; roll_pid];

    %% Dynamics/Plant
    % opts = odeset('MaxStep',1e-1);
    [t, x] = ode23(x_dot_fun, t0:dt:t0+t_inc, x0', [], u_c, sys, w_const);

    %% Storage
    if mode == 2
        x(:, 7:9) = wrapToPi(x(:, 7:9)); % wrapping Lambda between -pi to pi
    end
    x(:, 13:15) = wrapToPi(x(:, 13:15)); % wrapping alpha between 0 to 2pi
    thist = [thist; t(1:end-1)];
    xhist = [xhist; x(1:end-1, :)];
    u_chist = [u_chist, u_c];
    pids_hist = [pids_hist, pids];
    x0 = x(end, :)';
    t0 = t0 + t_inc;
    
    if ground(h) && t0 > 5
        disp(tspan(i));
        break;
    end
end

    %% Store Thrust and r w.r.t AOA
    toc;
    r_alpha_hist = [r_alpha_hist; r];
    T_alpha_hist = [T_alpha_hist; mean(u_chist(1:3,end-500:end),'all')];
    W_alpha_hist = [W_alpha_hist; mean(u_chist(4:6,end-500:end),'all')];
    
end

%% Store control metrics for performance evaluation 
u_chist_mode = [u_chist_mode; u_chist];

end

%% Performance Evaluation
for i = 1:numel(modes)
    u_chist = u_chist_mode(6*i-5:6*i,:);
    
    T_all = u_chist(1:3,:)./g;                  % Thrust over time (kgs)
    I_T = current_T(T_all);                     % Amperes (from experimental trend)
    avg_V_T = 17.5;                             % Volts
    % tot_P_T = sum(avg_V_T .* I_T, 1);           % Watts
    tot_P_T = sum(power_T(T_all), 1);           % Watts (from datasheet trend)
    
    tau_all =  u_chist(4:6,:);        % Torque on wings over time (Nm)
    I_tau = current_tau(abs(tau_all));          % Amperes
    avg_V_tau = 8;                              % Volts
    tot_P_tau = sum(avg_V_tau .* I_tau, 1);     % Watts
    
    tot_P = tot_P_T + tot_P_tau;                % Watts
    cum_P = cumsum(tot_P);                      % Watts
    cum_Wh = cum_P .* (t_inc/3600);             % Wh
    tot_flight_time = tot_Wh / mean(tot_P) * 60; % Minutes
    disp(tot_flight_time)

    % if i == 1
    %     save("var_tri", "t_span", "u_chist", "tot_flight_time", "cum_Wh")
    % elseif i == 2
    %     save("var_windmill", "t_span", "u_chist", "tot_flight_time", "cum_Wh")
    % elseif i == 3
    %     save("var_fixed", "t_span", "u_chist", "tot_flight_time", "cum_Wh")
    % end
end
return

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

for i = 1:18
    subplot(6,3,i)
    plot(thist, xplot(:, i))
    hold on
    if i == 1
        plot(tspan(1:length(deshist)), deshist(1, :))
    elseif i == 2
        plot(tspan(1:length(deshist)), deshist(2, :))
    elseif i == 3
        plot(tspan(1:length(deshist)), deshist(3, :))
    elseif i == 9 && mode == 1
        plot(tspan(1:length(deshist)), deshist(4, :))
    end
    xlabel("time (sec)");
    ylabel(x_name(i),'Interpreter','latex')
    grid on
end

%% Control Plots
set(groot,'defaultLineLineWidth',1)
set(groot,'defaultaxesFontSize',12)
fig2 = figure(2);
T_leg = ["T1", "T2", "T3"];
plot(tspan(1:length(u_chist)),u_chist(1:3,:)./g)
grid on
xlabel('time (sec)')
ylabel('kgs')
title('Thrust over time')
legend(T_leg)

fig3 = figure(3);
W_leg = ["W1", "W2", "W3"];
plot(tspan(1:length(u_chist)),u_chist(4:6,:))
grid on
xlabel('time (sec)')
ylabel('Nm')
title('Wing Torque over time')
legend(W_leg)

fig4 = figure(4);
plot(tspan(1:length(tot_P)), tot_P)
grid on
hold on
xlabel('time (sec)')
ylabel('Watts')
xline(8,"LineWidth",2,"Label","Tricopter","LabelVerticalAlignment","bottom")
xline(20,"LineWidth",2,"Label","Transition","LabelVerticalAlignment","bottom")
xline(26,"LineWidth",2,"Label","Windmill","LabelVerticalAlignment","bottom")
xl = [8,8,20,20,26,26,50,50];
yl = [0,700,700,0,0,700,700,0];
fill(xl(1:4),yl(1:4),'r','FaceAlpha',0.2)
fill(xl(3:6),yl(3:6),'b','FaceAlpha',0.2)
fill(xl(5:8),yl(5:8),'g','FaceAlpha',0.2)
xlim([8,tend])
ylim([0,700])
title('Power draw (Tri-copter to windmill mode transition)')

% fig5 = figure(5);
% pids_leg = ["h_pid", "pitch_pid", "roll_pid"];
% plot(tspan(1:length(u_chist)),pids_hist(2:3,:)')
% xlabel('time (sec)')
% ylabel('kgs')
% legend(pids_leg(2:3))

%% Performance Plots
mode_leg = ["Tri-copter", "Windmill", "Fixed-wing"];
for i = 1:3
    if i == 1
        load("var_tri", "tspan", "u_chist", "tot_flight_time", "cum_Wh")
    elseif i == 2
        load("var_windmill", "tspan", "u_chist", "tot_flight_time", "cum_Wh")
    elseif i == 3
        load("var_fixed", "tspan", "u_chist", "tot_flight_time", "cum_Wh")
    end

    fig6 = figure(6);
    plot(tspan(1:length(u_chist))/60, mean(u_chist(1:3,:)./g))
    grid on
    hold on
    xlabel('time (min)')
    ylabel('kgs')
    title('Thrust over time')
    legend(mode_leg)
    
    fig7 = figure(7);
    plot(tspan(1:length(u_chist))/60, mean(u_chist(4:6,:)))
    grid on
    hold on
    xlabel('time (min)')
    ylabel('Nm')
    title('Wing Torque over time')
    legend(mode_leg)
    
    fig8 = figure(8);
    plot(tspan(1:length(cum_Wh))/60, (cum_Wh/tot_Wh) * 100)
    grid on
    hold on
    xlabel('time (min)')
    ylabel('%')
    ylim([0,50])
    % yline(tot_Wh)
    title('Percent of Total energy consumed over time')
    legend(mode_leg)
end

%% Plot thrust vs AOA
fig9 = figure(9);
yyaxis("right")
plot(alphas,T_alpha_hist/g)
ylabel('kgs')
hold on
yline(0.52)
grid on
yyaxis("left")
plot(alphas, r_alpha_hist)
ylabel('rad/s')
xlabel('Wing Angle (deg)')
title('Thrust and Ang. Vel. w.r.t AOA')
legend('Ang. Vel.','Avg Thrust Windmill')

%% Animate
VIDEO = false;
video_title = "..\Videos\hover_drift";
fig10 = animate_system_seq(xhist, thist, VIDEO, video_title, t_inc);

%% Functions
function ground = ground(h)
ground = true;    
    if h>0
        ground = false;
    end
end

%% Notes
% Fixed wing mode isn't working; may be speed is not enough for lift