function Cm = Cmx(alpha)

sign_M = sign(alpha);
alpha = rad2deg(abs(alpha)); % change to deg from rad
if alpha >= 90
    alpha = 180 - alpha; % makes alpha negative
end

if alpha <= 18
    Cm = -0.0002 * alpha^3 + 0.0075 * alpha^2 - 0.088 * alpha + 0.1;
elseif alpha <= 47
    Cm = -0.0156 * alpha + 0.0385;
elseif alpha <= 90
    Cm = 9E-05 * alpha^2 - 0.0166 * alpha - 0.0361;
end

Cm = Cm * sign_M;
