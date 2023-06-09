function [kp, kd, ki] = run_gains(m, Jn, Je, Jd, Iw, max_T, vel_max, mode, alpha_des)

if mode == 1

    % Vertical Speed gains
    zeta_w = 0.707;
    a_w3 = 1/m;
    wn_w = sqrt(a_w3*max_T*sqrt(1-zeta_w^2)/vel_max);
    
    kp.w = wn_w^2/a_w3;
    kd.w = -1.1*(2*zeta_w*wn_w)/a_w3;
    ki.w = 0.3;
    
    % Height gains
    zeta_h = 0.707;
    a_h3 = 1/m;
    h_max = 0.8;
    wn_h = sqrt(a_h3*vel_max*sqrt(1-zeta_h^2)/h_max);
    
    kp.h = wn_h^2/a_h3;
    kd.h = -1.1*(2*zeta_h*wn_h)/a_h3;
    ki.h = 0.6;

    % Roll gains
    zeta_roll = 0.707;
    a_phi3 = 1/Jn;
    phi_max = deg2rad(10);
    tau_phi_max = 2;
    wn_roll = sqrt(a_phi3*tau_phi_max*sqrt(1-zeta_roll^2)/phi_max);
    
    kp.roll = wn_roll^2/a_phi3;
    kd.roll = 1.1*(2*zeta_roll*wn_roll)/a_phi3;
    ki.roll = 0;
    
    % Pitch gains
    zeta_pitch = 0.707;
    a_theta3 = 1/Je;
    theta_max = deg2rad(10);
    tau_theta_max = 2;
    wn_pitch = sqrt(a_theta3*tau_theta_max*sqrt(1-zeta_pitch^2)/theta_max);
    
    kp.pitch = wn_pitch^2/a_theta3;
    kd.pitch = 1.1*(2*zeta_pitch*wn_pitch)/a_theta3;
    ki.pitch = 0;

    % Yaw gains
    zeta_yaw = 0.707;
    a_psi3 = 1/Jd;
    yaw_max = deg2rad(30);
    tau_psi_max = 1;
    wn_yaw = sqrt(a_psi3*tau_psi_max*sqrt(1-zeta_yaw^2)/yaw_max);
    
    kp.yaw = wn_yaw^2/a_psi3;
    kd.yaw = 1.1*(2*zeta_pitch*wn_yaw)/a_psi3;
    ki.yaw = 0;
    
    % Alpha gains
    zeta_alpha = 0.707;
    a_alpha3 = 1/Iw;
    alpha_max = deg2rad(30);
    tau_alpha_max = 2;
    wn_alpha = sqrt(a_alpha3*tau_alpha_max*sqrt(1-zeta_alpha^2)/alpha_max);
    
    kp.alpha = wn_alpha^2/a_alpha3;
    kd.alpha = 1.1*(2*zeta_alpha*wn_alpha)/a_alpha3;
    ki.alpha = 0;
    
elseif mode == 2 || mode == 0
    
    % Vertical Speed gains
    zeta_w = 0.707;
    a_w3 = 1/m;
    w_max = 2 + ((2-4)/(45-75)) * (alpha_des-45);
    wn_w = sqrt(a_w3*max_T*sqrt(1-zeta_w^2)/w_max);
    
    kp.w = wn_w^2/a_w3;
    kd.w = -1.1*(2*zeta_w*wn_w)/a_w3;
    ki.w = 0.3;
    
    % Height gains
    zeta_h = 0.707;
    a_h3 = 1/m;
    h_max = 0.8 + ((0.8-1)/(45-75)) * (alpha_des-45);
    wn_h = sqrt(a_h3*w_max*sqrt(1-zeta_h^2)/h_max);
    
    kp.h = wn_h^2/a_h3;
    kd.h = -1.1*(2*zeta_h*wn_h)/a_h3;
    ki.h = 0.3;

    % Roll gains
    zeta_roll = 0.707;
    a_phi3 = 1/Jn;
    phi_max = deg2rad(30);
    tau_phi_max = 0.9;
    wn_roll = sqrt(a_phi3*tau_phi_max*sqrt(1-zeta_roll^2)/phi_max);
    
    kp.roll = wn_roll^2/a_phi3;
    kd.roll = 1.1*(2*zeta_roll*wn_roll)/a_phi3;
    ki.roll = 0;
    
    % Pitch gains
    zeta_pitch = 0.707;
    a_theta3 = 1/Je;
    theta_max = deg2rad(30);
    tau_theta_max = 1.15;
    wn_pitch = sqrt(a_theta3*tau_theta_max*sqrt(1-zeta_pitch^2)/theta_max);
    
    kp.pitch = wn_pitch^2/a_theta3;
    kd.pitch = 1.1*(2*zeta_pitch*wn_pitch)/a_theta3;
    ki.pitch = 0;
    
    % Yaw gains
    zeta_yaw = 0.707;
    a_psi3 = 1/Jd;
    yaw_max = deg2rad(30);
    tau_psi_max = 1;
    wn_yaw = sqrt(a_psi3*tau_psi_max*sqrt(1-zeta_yaw^2)/yaw_max);
    
    kp.yaw = wn_yaw^2/a_psi3;
    kd.yaw = 1.1*(2*zeta_pitch*wn_yaw)/a_psi3;
    ki.yaw = 0;

    % Alpha gains
    zeta_alpha = 0.707;
    a_alpha3 = 1/Iw;
    alpha_max = deg2rad(30);
    tau_alpha_max = 2;
    wn_alpha = sqrt(a_alpha3*tau_alpha_max*sqrt(1-zeta_alpha^2)/alpha_max);
    
    kp.alpha = wn_alpha^2/a_alpha3;
    kd.alpha = 1.1*(2*zeta_alpha*wn_alpha)/a_alpha3;
    ki.alpha = 0;

elseif mode == 3

    % Vertical Speed gains
    zeta_w = 0.707;
    a_w3 = 1/m;
    wn_w = sqrt(a_w3*max_T*sqrt(1-zeta_w^2)/vel_max);
    
    kp.w = wn_w^2/a_w3;
    kd.w = -1.1*(2*zeta_w*wn_w)/a_w3;
    ki.w = 0.3;
    
    % Height gains
    zeta_h = 0.707;
    a_h3 = 1/m;
    h_max = 0.8;
    wn_h = sqrt(a_h3*vel_max*sqrt(1-zeta_h^2)/h_max);
    
    kp.h = wn_h^2/a_h3;
    kd.h = -1.1*(2*zeta_h*wn_h)/a_h3;
    ki.h = 0.6;

    % Roll gains
    zeta_roll = 0.707;
    a_phi3 = 1/Jn;
    phi_max = deg2rad(10);
    tau_phi_max = 2;
    wn_roll = sqrt(a_phi3*tau_phi_max*sqrt(1-zeta_roll^2)/phi_max);
    
    kp.roll = wn_roll^2/a_phi3;
    kd.roll = 1.1*(2*zeta_roll*wn_roll)/a_phi3;
    ki.roll = 0;
    
    % Pitch gains
    zeta_pitch = 0.707;
    a_theta3 = 1/Je;
    theta_max = deg2rad(10);
    tau_theta_max = 2;
    wn_pitch = sqrt(a_theta3*tau_theta_max*sqrt(1-zeta_pitch^2)/theta_max);
    
    kp.pitch = wn_pitch^2/a_theta3;
    kd.pitch = 1.1*(2*zeta_pitch*wn_pitch)/a_theta3;
    ki.pitch = 0;
    
    % Yaw gains
    zeta_yaw = 0.707;
    a_psi3 = 1/Jd;
    yaw_max = deg2rad(30);
    tau_psi_max = 0.2;
    wn_yaw = sqrt(a_psi3*tau_psi_max*sqrt(1-zeta_yaw^2)/yaw_max);
    
    kp.yaw = wn_yaw^2/a_psi3;
    kd.yaw = 1.1*(2*zeta_pitch*wn_yaw)/a_psi3;
    ki.yaw = 0;

    % Alpha gains
    zeta_alpha = 0.707;
    a_alpha3 = 1/Iw;
    alpha_max = deg2rad(30);
    tau_alpha_max = 2;
    wn_alpha = sqrt(a_alpha3*tau_alpha_max*sqrt(1-zeta_alpha^2)/alpha_max);
    
    kp.alpha = wn_alpha^2/a_alpha3;
    kd.alpha = 1.1*(2*zeta_alpha*wn_alpha)/a_alpha3;
    ki.alpha = 0;
    
end

% Position Gains
kp.n = 0.075; kd.n = 0.15; ki.n = 0;
kp.e = 0.075; kd.e = 0.15; ki.e = 0;