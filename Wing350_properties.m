% Drone with 550mm Wing

%% Mass Properties
m = 1.528877;   % kg
Jn = 0.012001;  % kgm^2
Je = 0.012001;  % kgm^2
Jd = 0.022096;  % kgm^2

%% Wing Characteristics (From xflr5)
Iw = 0.000150;    % Inertia for Wing at rotation axis
l     = 0.29;         % Arm Length (m)
b     = l * 2;        % Wingspan (m)
AR    = 6.306;        % Aspect ratio
tr    = 1.429;         % Taper ratio
S     = (pi*l^2)/3; % 0.116/2;      % Wing Area (m^2)
y_MAC = (b / 6) * ((1 + 2 * tr) / (1 + tr)); % Wing location of MAC
MAC   = 0.086;        % Mean Aerodynamic Chord