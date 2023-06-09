function CL = CLx(alpha)

CL_alpha_0 = 0.0; % From xflr5
CL_alpha = 0.0879; % From xflr5

sign_L = sign(alpha);

alpha = rad2deg(abs(alpha)); % change to deg from rad
if alpha > 90
    alpha = 180 - alpha; % makes alpha negative
    sign_L = -sign_L;
end

if alpha <= 10
    CL = CL_alpha * alpha + CL_alpha_0;
elseif alpha <= 17
    CL = -0.0233 * alpha + 1.0186;
elseif alpha <= 56
    CL = -0.0008 * alpha^2 + 0.0714 * alpha - 0.438;
elseif alpha <= 90
    CL = -0.0258 * alpha + 2.4887;
end

CL = sign_L * CL;