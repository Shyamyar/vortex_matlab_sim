function LDR = LDR(alpha)

alpha = rad2deg(alpha); % change to deg from rad
LDR = CL(alpha) / CD(alpha);
% LDR = 0.001 * alpha^4 - 0.0031 * alpha^3 - ...
%        0.6185 * alpha^2 + 7.1064 * alpha - 5.7227; % Directly from xflr5
