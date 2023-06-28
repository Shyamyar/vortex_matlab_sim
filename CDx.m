function CD = CDx(alpha)
% CD = Cd0 + K * Cl^2;

alpha = rad2deg(abs(alpha)); % change to deg from rad
if alpha >= 90
    alpha = 180 - alpha; % makes alpha negative
end

if alpha < 15
    CD =0.0002888634*alpha^2 + 0.0012147336*alpha + 0.0092611732;
elseif alpha < 32
    CD =-0.0002401973*alpha^2 + 0.0287403309*alpha - 0.2519000443;
elseif alpha <= 90
    CD =-0.0001306885*alpha^2 + 0.0265938993*alpha  - 0.2967058959;
end
