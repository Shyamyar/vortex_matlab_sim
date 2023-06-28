function I = current_tau(tau)

% Hitec HS-7245MH
% Maximum Torque Range kg. / cm.	5.2 ~ 6.4
% Current Draw at Idle	            12 mA
% No Load Operating Current Draw	190 mA
% Stall Current Draw	            1,600 mA

max_tau = 6.4e-2 * 9.81;
max_I = 1600e-3;
I = (max_I/max_tau) .* tau + 12e-3;