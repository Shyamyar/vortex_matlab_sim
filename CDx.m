function CD = CDx(alpha)
% CD = Cd0 + K * Cl^2;

alpha = rad2deg(abs(alpha)); % change to deg from rad
if alpha >= 90
    alpha = 180 - alpha; % makes alpha negative
end

if alpha <= 15
    CD = 0.0008 * alpha^2 - 0.0031 * alpha + 0.0268;
elseif alpha <= 90
    CD = -0.0001 * alpha^2 + 0.0218 * alpha - 0.1354; % from xflr
end
