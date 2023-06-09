function Cm = Cm(alpha)

alpha = rad2deg(abs(alpha)); % change to deg from rad
if alpha <= 18
    Cm = -0.0002 * alpha^3 + 0.0075 * alpha^2 - 0.088 * alpha + 0.1;
elseif alpha <= 47
    Cm = -0.0156 * alpha + 0.0385;
else
    Cm = -6e-06 * alpha^3 + 0.0012 * alpha^2 - 0.0783 * alpha + 1;
end
