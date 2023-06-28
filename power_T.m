function P = power_T(T)

T = T .* 1000; % in grams
P = 0.0002 .* T.^2 + 0.1571 .* T - 1.3693;