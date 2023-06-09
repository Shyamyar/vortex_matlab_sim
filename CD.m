function CD = CD(alpha)
% CD = Cd0 + K * Cl^2;

alpha = rad2deg(abs(alpha)); % change to deg from rad
CD = -5e-6 * alpha^3 + 0.0004 * alpha^2 + 0.0048 * alpha; % from xflr
