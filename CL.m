function CL = CL(alpha)

CL_alpha_0 = 0.0; % From xflr5
CL_alpha = 0.0879; % From xflr5
alpha = rad2deg(alpha); % change to deg from rad
if alpha <= 10
    CL = CL_alpha * alpha + CL_alpha_0;
elseif alpha > 10 && alpha <= 17
    CL = -0.0233 * alpha + 1.0186;
else
    CL = -0.0008 * alpha^2 + 0.0714 * alpha - 0.438;
end
