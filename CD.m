function CD = CD(alpha)

alpha = rad2deg(alpha); % change to deg from rad
% CD = Cd0 + K * Cl^2;
CD = 0.0128 * alpha - 0.0124; % from xflr