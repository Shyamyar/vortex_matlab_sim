% Old Wing Wroking Simulation
%% Mass Properties 
m = 1.95;   % kg
Jn = 232225.1780 * 10^-6;  % kgm^2
Je = 232300.2777 * 10^-6;  % kgm^2
Jd = 454557.3979 * 10^-6;  % kgm^2

%% Wing Characteristics (From xflr5)
Iw = 645.99 * 10^-6;    % Inertia for Wing at rotation axis
l     = 0.47;         % Arm Length (m)
b     = l * 2;        % Wingspan (m)
AR    = 7.932;        % Aspect ratio
tr    = 1.37;         % Taper ratio
S     = (pi*l^2)/3; % 0.116/2;      % Wing Area (m^2)
y_MAC = (b / 6) * ((1 + 2 * tr) / (1 + tr)); % Wing location of MAC
MAC   = 0.124;        % Mean Aerodynamic Chord