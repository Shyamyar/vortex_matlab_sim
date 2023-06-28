function CL = CLx(alpha)

sign_L = sign(alpha);

alpha = rad2deg(abs(alpha)); % change to deg from rad
if alpha > 90
    alpha = 180 - alpha; % makes alpha negative
    sign_L = -sign_L;
end

if alpha < 15
    CL =-0.0025198237*alpha^2 + 0.1091744252*alpha - 0.0081683948;
elseif alpha < 32
    CL =0.0022324018*alpha^2 - 0.0965137835*alpha + 1.9029451107;
elseif alpha <= 90
    CL =-0.0007354709*alpha^2 + 0.0733573917*alpha - 0.5493489893;
end

CL = sign_L * CL;