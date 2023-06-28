function Cm = Cmx(alpha)

sign_M = sign(alpha);
alpha = rad2deg(abs(alpha)); % change to deg from rad
if alpha >= 90
    alpha = 180 - alpha; % makes alpha negative
end

if alpha < 15
    Cm =0.0012096563*alpha^2 - 0.0380895546*alpha + 0.0049156276;
elseif alpha < 32
    Cm =-0.0006131512*alpha^2 + 0.0165526915*alpha - 0.4062500492;
elseif alpha <= 90
    Cm =0.000175835*alpha^2 - 0.0286469325*alpha + 0.2451278964;
end

Cm = Cm * sign_M;
